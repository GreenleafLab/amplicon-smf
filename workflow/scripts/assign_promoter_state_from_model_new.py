# input a series of single molecule matrices and a pretrained model
# for now only supports kmeans, we might be able to generalize this in the future?
# returns a dataframe with index = molecule ID and cluster assignment and sample and amplicon info
# also returns a dataframe with just the promoter GpCs in the same order if necessary
# optionally in the future can add details about the clustering and also additional metadata columns? prob want to do this after the fact though
# optionally can take in an offset to account for differences in sequencing primers between the data and the model
# note that I don't actually think we need to do this, since as long as we align to the same amplicon the coordinates of the promoter will be the same
# which works for us now given that we revcomp the sequences, but maybe not in the future?
# for now only works with using a single basedir and a list of samples and amplicons, in the future can customize, but this basically works for the current output of the pipeline

import pandas as pd
import numpy as np
import os
from os import path
from matplotlib import pyplot as plt
import argparse
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import seaborn as sns
from sklearn.cluster import KMeans
import pickle
from copy import deepcopy
import gc

from common import load_single_molecule_matrix, adjust_gcgs, annotate_40_clusters, assign_cluster, assign_promoter_footprints

cluster_order = [23,16,1,6,36,34,12,4,29,7,22,18,35,19,37,32,3,9,13,14,26,15,2,10,25,31,0,17,28,24,20,11,33,21,27,39,5,8,38,30]
new_cluster_order = [1,4,5,13,0,6,16,11,17,18,2,3,15,12,8,9,19,14,10,7]
promoter_on = [24,11,20,33,21,39,27,5,8,38,30]

tata = (155,170)
tss = (120,135)
pause = (80,100)

promoter_regions = {
    'tata': tata,
    'tss': tss,
    'pause': pause,
}

def get_all_protection_streaks(methyls, positions, length_method='first_to_last'):
    '''
    This function finds all streaks of protection in the molecule.
    
    Parameters:
    -----------
    methyls : list
        List of methylation values (0 or 1)
    positions : list
        List of genomic positions corresponding to methyls
    length_method : str
        Method for calculating streak length:
        - 'first_to_last': Distance from first 1 to last 1 (original)
        - 'interpolate': Midpoint between flanking 0s and 1s
        - 'flanking_zeros': Distance from last 0 before to first 0 after
    
    Returns:
    --------
    list of tuples: (length, start, end) sorted by length descending
    '''
    
    # Select the length calculation function once
    if length_method == 'first_to_last':
        calc_length = lambda start, end, last_zero, next_zero: (end - start, start, end)
    elif length_method == 'interpolate':
        calc_length = lambda start, end, last_zero, next_zero: (
            (end + next_zero) / 2 - (last_zero + start) / 2,
            (last_zero + start) / 2,
            (end + next_zero) / 2
        )
    elif length_method == 'flanking_zeros':
        calc_length = lambda start, end, last_zero, next_zero: (
            next_zero - last_zero,
            last_zero,
            next_zero
        )
    else:
        raise ValueError(f"Unknown length_method: {length_method}. "
                        f"Use 'first_to_last', 'interpolate', or 'flanking_zeros'")
    
    methyls = methyls + [0]
    positions = positions + [np.inf]
    
    streak = False
    start_pos = positions[0] - 1
    end_pos = positions[0] - 1
    last_zero_pos = positions[0] - 1
    streaks = []

    for idx in range(len(methyls)):
        if streak:
            if methyls[idx] == 1:
                end_pos = positions[idx]
            else:
                # End of streak - use the selected calculation function
                streak = False
                length, boundary_start, boundary_end = calc_length(
                    start_pos, end_pos, last_zero_pos, positions[idx]
                )
                streaks.append((length, boundary_start, boundary_end))
                last_zero_pos = positions[idx]
        else:
            if methyls[idx] == 1:
                start_pos = positions[idx]
                end_pos = positions[idx]
                streak = True
            else:
                last_zero_pos = positions[idx]
            
    return sorted(streaks, key=lambda x: x[0], reverse=True)

# # change this back
# def get_all_protection_streaks(methyls, positions):
#     '''
#     This function basically finds all streaks of protection in the molecule
#     Note it doesn't take into account streaks that are interrupted by a random methylation (ugh..)
#     '''

#     # append a 0 to methyls so we always finish with a nonmethyl (so we don't need edge cases)
#     methyls = methyls + [0]
#     # also need to do likewise to positions
#     positions = positions + [np.inf]
#     streak = False
#     start = -1
#     end = -1
#     length = 0
    
#     streaks = []

#     for idx in range(len(methyls)):
#         if streak:
#             if methyls[idx] == 1:
#                 end = positions[idx]
#                 length = end - start
#             else:
#                 streak = False
#                 streaks.append((length,start,end))
#                 start = positions[idx]
#                 end = positions[idx]
#                 length = 0
#         else:
#             if methyls[idx] == 1:
#                 start = positions[idx]
#                 end = positions[idx]
#                 length = 0
#                 streak = True
#             else:
#                 pass
            
#     return sorted(streaks, key=lambda x: x[0], reverse=True)

def remove_nucs_from_promoter_matrix(mat, promoter_region_threshold_high, nuc_length_threshold, 
                                     min_3prime_nuc_length, tss, promoter_low, promoter_high, 
                                     length_method='interpolate'):
    '''
    This function finds all nucleosomes and removes them from the promoter matrix.
    
    Returns:
    --------
    no_nuc_mat : DataFrame
        Matrix with nucleosomes removed, cropped to promoter region
    nuc_over_tss_list : list
        Boolean - does a nuc overlap TSS
    nuc_over_promoter_list : list
        Boolean - does a nuc overlap promoter region
    has_plus1_list : list
        Boolean - is there a nucleosome downstream of TSS (+1 position)
    has_minus1_list : list
        Boolean - is there a nucleosome upstream of TSS (-1 position)
    plus1_distance_list : list
        Distance from TSS to upstream edge of first downstream nuc (NaN if none)
    '''

    positions = mat.columns.tolist()
    no_nuc_mat = mat.copy()
    
    # figure out the phantom bases weirdness
    phantom_length = nuc_length_threshold - min_3prime_nuc_length
    molecule_end = min(mat.columns.tolist())
    ghost_start = molecule_end - phantom_length
    fake_bases = list(range(ghost_start, molecule_end, 10))
    fake_protections = [1 for _ in fake_bases]

    nuc_over_tss_list = []
    nuc_over_promoter_list = []
    has_plus1_list = []
    has_minus1_list = []
    plus1_distance_list = []
    
    # iterate through mat
    for idx in mat.index:
        nuc_over_tss = False
        nuc_over_promoter = False
        plus1_distance = np.nan
        # has_minus1 = False

        single_molecule_signal = mat.loc[idx]

        this_methyl_positions = fake_bases + single_molecule_signal.index.tolist()
        methyls = fake_protections + single_molecule_signal.tolist()

        protection_streaks = get_all_protection_streaks(methyls, this_methyl_positions, 
                                                        length_method=length_method)

        # iterate through all nucleosomes
        for streak in protection_streaks:
            length = streak[0]
            lower_bound = streak[1]  # Smaller coordinate (downstream)
            upper_bound = streak[2]  # Larger coordinate (upstream)

            if length >= nuc_length_threshold:
                # it's a nuc - remove it from matrix
                cols = [c for c in positions if c >= lower_bound and c <= upper_bound]
                no_nuc_mat.loc[idx, cols] = 0

                # Check if nuc overlaps TSS
                if tss > lower_bound and tss < upper_bound:
                    nuc_over_tss = True
                
                # Check if nuc overlaps promoter
                if lower_bound < promoter_high and upper_bound > promoter_low:
                    nuc_over_promoter = True
                
                # Find +1 nucleosome (downstream of TSS = lower coordinate)
                if upper_bound <= tss:
                    # Nucleosome is completely downstream of TSS
                    distance = tss - upper_bound
                    if np.isnan(plus1_distance) or distance < plus1_distance:
                        plus1_distance = distance
                
                # # Check for -1 nucleosome (upstream of TSS = higher coordinate)
                # if lower_bound >= tss:
                #     # Nucleosome is completely upstream of TSS
                #     has_minus1 = True

        # Determine if +1 present
        has_plus1 = not np.isnan(plus1_distance)
        
        nuc_over_tss_list.append(nuc_over_tss)
        nuc_over_promoter_list.append(nuc_over_promoter)
        has_plus1_list.append(has_plus1)
        # has_minus1_list.append(has_minus1)
        plus1_distance_list.append(plus1_distance)
    
    return (no_nuc_mat[[c for c in no_nuc_mat.columns if c <= promoter_region_threshold_high]], 
            nuc_over_tss_list, 
            nuc_over_promoter_list,
            has_plus1_list,
            # has_minus1_list,
            plus1_distance_list)

# def remove_nucs_from_promoter_matrix(mat, promoter_region_threshold_high, nuc_length_threshold, min_3prime_nuc_length, tss, promoter_low, promoter_high, length_method='interpolate'):
#     '''
#     This function basically first finds all nucleosomes (SLOW but who cares) and then removes them from the promoter matrix
#     So that we can cluster on only the underlying mol bio at the promoter
#     And then we can also use this to just do a quick check of whether there is a nuc that overlaps the promoter (or TSS)
#     '''

#     positions = mat.columns.tolist()
#     no_nuc_mat = mat.copy()
    
#     # figure out the phantom bases weirdness
#     phantom_length = nuc_length_threshold - min_3prime_nuc_length
#     molecule_end = min(mat.columns.tolist())
#     ghost_start = molecule_end - phantom_length
#     fake_bases = list(range(ghost_start, molecule_end, 10))
#     fake_protections = [1 for _ in fake_bases]

#     nuc_over_tss_list = []
#     nuc_over_promoter_list = []
    
#     # iterate through mat
#     for idx in mat.index:
#         nuc_over_tss, nuc_over_promoter = False, False
#         plus1_nuc, minus1_nuc = False, False
#         dist_to_plus1 = np.nan

#         single_molecule_signal = mat.loc[idx]

#         this_methyl_positions = fake_bases + single_molecule_signal.index.tolist()
#         methyls = fake_protections + single_molecule_signal.tolist()

#         protection_streaks = get_all_protection_streaks(methyls, this_methyl_positions, length_method=length_method)

#         # iterate through all nucleosomes
#         for streak in protection_streaks:
#             length = streak[0]
#             lower_bound = streak[1]
#             upper_bound = streak[2]

#             if length >= nuc_length_threshold:
#                 # it's a nuc 
#                 cols = [c for c in positions if c >= lower_bound and c <= upper_bound]
#                 no_nuc_mat.loc[idx,cols] = 0

#                 # now, add other annotations
#                 if tss > lower_bound and tss < upper_bound: # nuc over tss
#                     # nuc over tss = true
#                     nuc_over_tss = True
#                 if lower_bound < promoter_high and upper_bound > promoter_low: # nuc over promoter
#                     # nuc over promoter = true
#                     nuc_over_promoter = True
#                 # # now ask about +1, and -1 nucs
#                 # if upper_bound < tss: # +1 nuc
#                 #     # make sure to check directionality/signs here, I probably messed this up tbh
#                 #     plus1_nuc = True
#                 #     # now figure out the distance from the tss to the +1 if it exists (this is nested, already know it's +1)

#         nuc_over_tss_list.append(nuc_over_tss)
#         nuc_over_promoter_list.append(nuc_over_promoter)
    
#     return no_nuc_mat[[c for c in no_nuc_mat.columns if c <= promoter_region_threshold_high]], nuc_over_tss_list, nuc_over_promoter_list

# def plot_promoters(mat, no_nuc_mat, plots_file, title):
#     fig, axs = plt.subplots(1, 2)
#     axs = axs.ravel()

#     mapping = {1: 0.0, 0: 0.7, -1: 1.0 - np.finfo(np.float32).eps}

#     A = mat.replace(mapping).to_numpy(dtype=np.float32, copy=False)
#     B = no_nuc_mat.replace(mapping).to_numpy(dtype=np.float32, copy=False)

#     axs[0].imshow(A, aspect='auto', cmap='gray', interpolation='none', vmin=0.0, vmax=1.0)
#     axs[1].imshow(B, aspect='auto', cmap='gray', interpolation='none', vmin=0.0, vmax=1.0)

#     axs[0].invert_xaxis()
#     axs[1].invert_xaxis()

#     axs[0].set_title('Original')
#     axs[1].set_title('Nucs Removed')
    
#     # Figure-level title
#     fig.suptitle(title, fontsize=14, fontweight='bold')
    
#     # Adjust layout to prevent title overlap
#     fig.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space at top for suptitle

#     plots_file.savefig(fig)
#     plt.close(fig)
#     del A, B, fig, axs

def assign_molecular_states(annot_df, state_map_dict):
    states_to_return = pd.DataFrame(index=annot_df.index)

    # do operations on tata/tbp/paused etc
    states_to_return['nuc'] = annot_df['nuc_over_promoter']
    states_to_return['promoter_nuc_free'] = ~annot_df['nuc_over_promoter']
    states_to_return['pic'] = annot_df['tata'] & annot_df['tss'] & ~annot_df['pause'] & states_to_return['promoter_nuc_free']
    states_to_return['tbp'] = annot_df['tata'] & ~annot_df['tss'] & ~annot_df['pause'] & states_to_return['promoter_nuc_free']
    states_to_return['only_pause'] = annot_df['pause'] & ~annot_df['tata'] & states_to_return['promoter_nuc_free']

    # do operations on the new cluster column
    states_to_return['cluster_category'] = annot_df['new_cluster'].map(state_map_dict)
    states_to_return['state'] = states_to_return['cluster_category']
    states_to_return.loc[states_to_return['promoter_nuc_free'] == False,'state'] = 'nuc'

    return states_to_return

def plot_promoters(mat, no_nuc_mat, plots_file, title, max_rows=10000):
    original_n = len(mat)
    
    # Subsample if needed
    if len(mat) > max_rows:
        print(f"  Subsampling {len(mat)} -> {max_rows} for plotting")
        sample_idx = np.random.choice(mat.index, size=max_rows, replace=False)
        sample_idx_sorted = [idx for idx in mat.index if idx in sample_idx]
        mat = mat.loc[sample_idx_sorted]
        no_nuc_mat = no_nuc_mat.loc[sample_idx_sorted]
    
    fig, axs = plt.subplots(1, 2, figsize=(12, 8))
    axs = axs.ravel()

    mapping = {1: 0.0, 0: 0.7, -1: 1.0 - np.finfo(np.float32).eps}

    A = mat.replace(mapping).to_numpy(dtype=np.float32, copy=True)
    B = no_nuc_mat.replace(mapping).to_numpy(dtype=np.float32, copy=True)

    axs[0].imshow(A, aspect='auto', cmap='gray', interpolation='none', vmin=0.0, vmax=1.0)
    axs[1].imshow(B, aspect='auto', cmap='gray', interpolation='none', vmin=0.0, vmax=1.0)

    axs[0].invert_xaxis()
    axs[1].invert_xaxis()
    
    axs[0].set_title('Original')
    axs[1].set_title('Nucs Removed')
    
    # Show if subsampled
    title_text = f"{title} (n={original_n})"
    if original_n > max_rows:
        title_text += f" [showing {max_rows}]"
    fig.suptitle(title_text, fontsize=14, fontweight='bold')
    
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    plots_file.savefig(fig, dpi=75)
    plt.close(fig)
    
    del A, B, fig, axs

def compare_promoter_accessibility(annot_df, plot_dir, prefix):
    fig, axs = plt.subplots(1,2, figsize=(10, 5), constrained_layout=True)
    axs = axs.ravel()

    to_plot = ( 
        annot_df.groupby(['sample','amplicon'])[['promoter_cluster_on','nuc_over_tss','nuc_over_promoter']]
        .mean()
        .reset_index() 
    ) 

    for ax, ycol in zip(axs, ['nuc_over_tss', 'nuc_over_promoter']):
        x = to_plot['promoter_cluster_on'].to_numpy()
        y = to_plot[ycol].to_numpy()

        hb = ax.hexbin(x, y, gridsize=60, bins='log', mincnt=1, cmap='magma')
        # fig.colorbar(hb, ax=ax, label='count (log)')

        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_box_aspect(1)
        ax.set_xlabel('promoter_cluster_on')
        ax.set_ylabel(ycol)
        ax.grid(True, alpha=0.15)

    fig.savefig(path.join(plot_dir, f'{prefix}.promoter_accessibility_calling.pdf'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assigns promoter state to a single-molecule matrix given a model')
    parser.add_argument("--input_basedir", dest="basedir", type=str, help="Path to highlevel folder containing multiple single molecule matrices")
    parser.add_argument("--samples", type=str, nargs='+', help="Space-separated list of sample names")
    parser.add_argument("--amplicons", type=str, nargs='+', help="Space-separated list of amplicon names")
    parser.add_argument("--output", type=str, help="Output prefix")
    parser.add_argument("--output_dir", type=str, default='.', help="Path to output folder")
    parser.add_argument("--plot_dir", type=str, default=None, help="Path to plot folder")
    parser.add_argument("--model", type=str, default="/oak/stanford/groups/wjg/bgrd/scripts_share/220623_kmeans_40clusters_promotermodel.pkl", help="Path to .pkl k-means model file to use")
    parser.add_argument("--new_model", type=str, default="/oak/stanford/groups/wjg/bgrd/papers/ad_smf/data/samplesheets/2602223_new_model_v2.pkl", help="Path to .pkl k-means model file to use")
    parser.add_argument("--cluster_to_state_dict", type=str, default="/oak/stanford/groups/wjg/bgrd/papers/ad_smf/data/samplesheets/260223_promoter_state_annotations_v1.txt", help="Path to txt file with annotations for new_cluster")
    parser.add_argument("--promoter_region_bound_low", type=int, default=0, help="Lower bound of promoter region (to subset methyl positions)")
    parser.add_argument("--promoter_region_bound_high", type=int, default=240, help="Upper bound of promoter region (to subset methyl positions)")
    parser.add_argument("--no_dedup", action='store_true', help="Whether to use the deduped input or not (defaults to dedup)")
    parser.add_argument("--dont_correct", action='store_true', help="Whether to do the CpG correction in the promoter or not (defaults to correcting)")
    parser.add_argument("--min_nuc_len", type=int, default=100, help="How long a stretch of protection has to be for us to treat it as a nucleosome")
    parser.add_argument("--min_3p_overlap", type=int, default=18, help="How many bases at the end of the molecule have to be protected to treat it as a nucleosome")
    parser.add_argument("--promoter_bound_low", type=int, default=136, help="Lower bound of minCMV (to check for nucs) = TSS")
    parser.add_argument("--promoter_bound_high", type=int, default=195, help="Upper bound of minCMV (to check for nucs)")
    parser.add_argument("--length_method", type=str, default='interpolate', choices=['first_to_last', 'interpolate', 'flanking_zeros'], help="Method for calculating protection streak length")

    args = parser.parse_args()

    promoter_dfs = []
    promoter_dfs_no_nuc = []
    annotation_dfs = []

    # iterate through all (sample, amplicon) tuples
    for sample in args.samples:
        print(sample)
        for amplicon in args.amplicons:
            # read in matrix
            p = path.join(args.basedir, sample, 'matrices', '{}.{}.full_unclustered.matrix'.format(sample, amplicon)) if args.no_dedup else path.join(args.basedir, sample, 'matrices', '{}.{}.dedup.full_unclustered.matrix'.format(sample, amplicon))
            if path.isfile(p):
                # load matrix 
                mat = load_single_molecule_matrix(p, every_other=args.no_dedup)

                # reset index with additional info and cast to smaller dtype
                mat.index = ['{}-{}-{}'.format(sample, amplicon, i) for i in mat.index]
                mat = mat.astype('int8', copy=False)

                # correct promoter positions that can be ambiguous
                if not args.dont_correct:
                    mat = adjust_gcgs(mat, 54, 50) # 4/3 -- note that for the new clustering with the every_other=False the numbers are off by 1
                    mat = adjust_gcgs(mat, 130, 123) # 9/8

                # extract nucleosomes from the matrix
                # this also crops to just promoter GpCs and returns whether there is a nuc over the promoter/tss
                no_nuc_mat, nuc_over_tss, nuc_over_promoter, has_plus1, plus1_dist = remove_nucs_from_promoter_matrix(
                    mat, 
                    args.promoter_region_bound_high, 
                    args.min_nuc_len, 
                    args.min_3p_overlap, 
                    args.promoter_bound_low, 
                    args.promoter_bound_low, 
                    args.promoter_region_bound_high,
                    length_method=args.length_method
                )

                # no_nuc_mat, nuc_over_tss, nuc_over_promoter = remove_nucs_from_promoter_matrix(mat, args.promoter_region_bound_high, args.min_nuc_len, args.min_3p_overlap, args.promoter_bound_low, args.promoter_bound_low, args.promoter_region_bound_high, length_method=args.length_method)
                
                # store the results 
                mat = mat[[c for c in mat.columns if c <= args.promoter_region_bound_high]]

                promoter_dfs.append(mat)
                promoter_dfs_no_nuc.append(no_nuc_mat)

                # instantiate annotation df with metadata
                annot_df = pd.DataFrame(index=mat.index)

                annot_df['sample'] = sample
                annot_df['amplicon'] = amplicon

                annot_df['nuc_over_tss'] = nuc_over_tss
                annot_df['nuc_over_promoter'] = nuc_over_promoter
                annot_df['has_plus1_nuc'] = has_plus1
                # annot_df['has_minus1_nuc'] = has_minus1
                annot_df['plus1_distance'] = plus1_dist

                annotation_dfs.append(annot_df)
            else:
                continue
                # print(sample, amplicon)
                # continue # print(p)

    # join all the sub tables
    print('start building tables')
    all_promoter_matrix = pd.concat(promoter_dfs)
    all_promoter_matrix_no_nuc = pd.concat(promoter_dfs_no_nuc)
    all_annotation_matrix = pd.concat(annotation_dfs)

    all_annotation_matrix['sample'] = all_annotation_matrix['sample'].astype('category')
    all_annotation_matrix['amplicon'] = all_annotation_matrix['amplicon'].astype('category')
    print('end building tables')

    print(len(all_annotation_matrix))
    all_annotation_matrix.info(memory_usage='deep')

    count = np.sum(all_promoter_matrix.to_numpy() == -1)
    print(f"Total -1 values: {count}")

    # annotate promoters (old way) 
    print('start cluster')
    clusters = assign_cluster(all_promoter_matrix, args.model)
    all_annotation_matrix['cluster'] = pd.Categorical(clusters, categories=cluster_order, ordered=True)
    print('end cluster')

    # annotate the clusters 
    print('start annotate')
    promoter_is_on = annotate_40_clusters(all_annotation_matrix, promoter_on)
    all_annotation_matrix['promoter_cluster_on'] = promoter_is_on
    print('end annotate')

    # add kmeans from new way
    print('start new kmeans')
    new_clusters = assign_cluster(all_promoter_matrix_no_nuc, args.new_model)
    all_annotation_matrix['new_cluster'] = pd.Categorical(new_clusters, categories=new_cluster_order, ordered=True)
    print('end new kmeans')

    # classify based on pre-determined footprints
    print('start classify footprints')
    promoter_states_df = assign_promoter_footprints(all_promoter_matrix_no_nuc, promoter_regions)
    all_annotation_matrix = all_annotation_matrix.join(promoter_states_df)
    print('end classify footprints')

    # now assign states like pic/pause etc
    # note -- for now we are going to do this on the called footprint positions, and some logic around that
    # but theoretically we should also use the clustering results from the k means 
    print('start annotate actual states')
    cluster_to_state_dict = pd.read_table(args.cluster_to_state_dict).set_index('cluster')
    cluster_to_state_dict = cluster_to_state_dict['cluster_category'].to_dict()
    state_annot_df = assign_molecular_states(all_annotation_matrix, cluster_to_state_dict)
    all_annotation_matrix = all_annotation_matrix.join(state_annot_df)
    print('end annotate actual states')

    # write to file
    print('begin write')
    all_promoter_matrix.to_csv(path.join(args.output_dir,'{}_promoter_methylation_data.txt.gz'.format(args.output)), sep='\t', index=True, header=True)
    all_promoter_matrix_no_nuc.to_csv(path.join(args.output_dir,'{}_promoter_methylation_data_nucs_removed.txt.gz'.format(args.output)), sep='\t', index=True, header=True)
    all_annotation_matrix.to_csv(path.join(args.output_dir,'{}_promoter_annotations.txt.gz'.format(args.output)), sep='\t', index=True, header=True)
    print('end write')

    # also now write the per amplicon aggregated stats for use later
    # need to basically define the states i'm gonna use here for promoters later and then get a percentage for each
    print('start amplicon level')
    def generate_clust_string(series):
        """
        Custom aggregator to create the cluster proportions string.
        """
        counts = series.value_counts(normalize=True).sort_index()
        return ','.join([f"{c}:{freq:.4f}" for c, freq in counts.items()])

    per_amp_df = all_annotation_matrix.groupby(['sample', 'amplicon']).agg(
        frac_promoter_open=('promoter_nuc_free', 'mean'),
        molecules=('cluster', 'count'),
        tbp=('tbp', 'mean'),
        pic=('pic', 'mean'),
        pause=('pause', 'mean'),
        cluster_proportions=('state', generate_clust_string)
    ).reset_index()

    # promoter_annotations_groupby = promoter_annotations.groupby(['sample','amplicon'])[['promoter_nuc_free','tbp','pic','pause','only_pause']].mean().reset_index()

    # promoter_annotations_groupby['fraction_tbp'] = promoter_annotations_groupby['tbp'] / promoter_annotations_groupby['promoter_nuc_free']
    # promoter_annotations_groupby['fraction_pic'] = promoter_annotations_groupby['pic'] / promoter_annotations_groupby['promoter_nuc_free']
    # promoter_annotations_groupby['fraction_pause'] = promoter_annotations_groupby['pause'] / promoter_annotations_groupby['promoter_nuc_free']
    # promoter_annotations_groupby['fraction_only_pause'] = promoter_annotations_groupby['only_pause'] / promoter_annotations_groupby['promoter_nuc_free']
    # promoter_annotations_groupby['pic_over_tbp'] = promoter_annotations_groupby['pic'] / promoter_annotations_groupby['tbp']
    # promoter_annotations_groupby['paused_over_pic'] = promoter_annotations_groupby['pause'] / promoter_annotations_groupby['pic']

    per_amp_df.to_csv(path.join(args.output_dir,'{}_promoter_annotations_amplicon_level.txt.gz'.format(args.output)), sep='\t', index=False, header=True)
    print('end amplicon level')

    # also do sample level, which is just 6x frac open!
    print('start samp level')
    amp = 'opJS4_6x_TetO_21bp_no_CG'
    per_samp_df = per_amp_df.loc[per_amp_df.amplicon==amp,['sample', 'frac_promoter_open']]
    per_samp_df.rename(columns={'frac_promoter_open': 'frac_promoter_open_6x'}).to_csv(path.join(args.output_dir,'{}_promoter_annotations_sample_level.txt.gz'.format(args.output)), sep='\t', index=False, header=True)
    print('end samp level')

    print('start plot')
    # optionally plot
    if args.plot_dir:
        # unique_samples = all_annotation_matrix['sample'].unique()
        # print(f"Plotting {len(unique_samples)} samples")
        for samp in args.samples:
            print(samp)

            df1 = all_annotation_matrix[all_annotation_matrix['sample'] == samp]
            pdf_path = path.join(args.plot_dir, f'{samp}.promoter_plots.pdf')

            with PdfPages(pdf_path) as plots:
                # unique_amplicons = df1['amplicon'].unique()

                for amp in args.amplicons:
                    # print(amp)
                    df2 = df1[df1['amplicon'] == amp]
                    idxs = df2.sort_values('cluster').index

                    # print(len(df2))

                    plot_promoters(all_promoter_matrix.loc[idxs], all_promoter_matrix_no_nuc.loc[idxs], plots, amp)

                    # Force cleanup after each plot
                    plt.close('all')
                    del df2, idxs
                    gc.collect()

            # Cleanup after each sample
            del df1
            gc.collect()

        # also compare the different ways of calling promoter accessibility
        compare_promoter_accessibility(all_annotation_matrix, args.plot_dir, args.output)

    print('end plot')

