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

from common import load_single_molecule_matrix


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assigns promoter state to a single-molecule matrix given a model')
    parser.add_argument("--input_basedir", dest="basedir", type=str, help="Path to highlevel folder containing multiple single molecule matrices")
    parser.add_argument("--samples", dest="samples", type=str, nargs='+', help="Space-separated list of sample names")
    parser.add_argument("--amplicons", dest="amplicons", type=str, nargs='+', help="Space-separated list of amplicon names")
    parser.add_argument("--output", dest="output", type=str, help="Output prefix")
    parser.add_argument("--output_dir", dest="output_dir", type=str, default='.', help="Path to output folder")
    parser.add_argument("--model", dest="model", type=str, default="/oak/stanford/groups/wjg/bgrd/scripts_share/220623_kmeans_40clusters_promotermodel.pkl", help="Path to .pkl k-means model file to use")
    parser.add_argument("--promoter_bound_low", dest="promoter_bound_low", type=int, default=0, help="Lower bound of promoter region (to subset methyl positions)")
    parser.add_argument("--promoter_bound_high", dest="promoter_bound_high", type=int, default=240, help="Upper bound of promoter region (to subset methyl positions)")

    args = parser.parse_args()

    promoter_dfs = []
    annotation_dfs = []

    # iterate through all (sample, amplicon) tuples
    for sample in args.samples:
        for amplicon in args.amplicons:
            # read in matrix
            p = path.join(args.basedir, sample, 'matrices', '{}.{}.full_unclustered.matrix'.format(sample, amplicon))
            if path.isfile(p):
                print(p)
                mat = load_single_molecule_matrix(p)

                # reset index with additional info
                mat.index = ['{}-{}-{}'.format(sample, amplicon, i) for i in mat.index]

                # extract just the promoter GpCs
                promoter_mat = mat[[c for c in mat.columns if c > args.promoter_bound_low and c < args.promoter_bound_high]]
                promoter_dfs.append(promoter_mat.copy())

                # instantiate annotation df with metadata
                annot_df = pd.DataFrame(index=mat.index)
                annot_df['sample'] = sample
                annot_df['amplicon'] = amplicon
                annotation_dfs.append(annot_df.copy())
            else:
                print(p)

    # join all the sub tables
    all_promoter_matrix = pd.concat(promoter_dfs)
    all_annotation_matrix = pd.concat(annotation_dfs)

    # load the model 
    with open(args.model, "rb") as m:
        kmeans_model = pickle.load(m)

    # cluster all of the promoters together
    promoter_clusters = kmeans_model.predict(all_promoter_matrix)
    all_annotation_matrix['cluster'] = promoter_clusters

    # write to file
    all_promoter_matrix.to_csv(path.join(args.output_dir,'{}_promoter_methylation_data.txt'.format(args.output)), sep='\t', index=True, header=True)
    all_annotation_matrix.to_csv(path.join(args.output_dir,'{}_promoter_annotations.txt'.format(args.output)), sep='\t', index=True, header=True)

