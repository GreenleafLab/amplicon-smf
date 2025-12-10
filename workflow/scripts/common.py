import pandas as pd
import numpy as np
# import os
# from os import path
from matplotlib import pyplot as plt
# from Bio import SeqIO
# import argparse
from matplotlib.backends.backend_pdf import PdfPages
plt.switch_backend('agg')
import seaborn as sns
# from Bio import SeqIO
# import re
from matplotlib.patches import Circle, Ellipse
from scipy.cluster.hierarchy import fcluster, fclusterdata, linkage, dendrogram 

# Note: assumes the input plot is a PdfPages (and so calls the savefig method on it)
# we can probably make this generic for other plotting types later (i.e. if pdfpages: plots.savefig() vs. ...)
def plot_bulk_smf_trace(df, plots, title='', loc_col='start', smf_col='SMF', ylim=None, ylabel='%SMF', fillbetween=None):
    '''
    Helper function for plotting a SMF trace across a locus
    Input is a dataframe with contig, position, and SMF%
    Returns a plot to the PdfPages object
    '''
    fig, ax = plt.subplots()

    df.plot(loc_col, smf_col, ls='-', ax=ax, c='k', ms=10, marker='.', legend=None, title=title)

    if ylim:
        ax.set_ylim(ylim)
    ax.set_xlabel('Amplicon position (bp)')
    ax.set_ylabel(ylabel)

    if fillbetween:
        add_fillbetween(ax, fillbetween)
        plt.legend()

    plt.tight_layout()
    plots.savefig()
    plt.close()

def plot_bulk_smf_trace_from_matrix(ser, plots, title='', ylim=None, ylabel='%SMF', fillbetween=None):
    '''
    Helper function for plotting a SMF trace across a locus
    Input is a series computed by taking the average down the columns from a single molecule matrix
    Returns a plot to the PdfPages object
    '''
    fig, ax = plt.subplots()

    ser.plot(ls='-', marker='.', ms=10, c='k', legend=None, title=title)

    if ylim:
        ax.set_ylim(ylim)
    ax.set_xlabel('Amplicon position (bp)')
    ax.set_ylabel(ylabel)

    if fillbetween:
        add_fillbetween(ax, fillbetween)
        plt.legend()

    plt.tight_layout()
    plots.savefig()
    plt.close()

def load_tfbs_positions(f):
    '''
    Loads a file of positions to color on the bulk SMF plots
    File looks like:
    >amplicon
    101,105,name
    299,420,name
    ...
    >amplicon2
    ...
    Returns dict of amplicon->list of (start,end) tuples
    '''
    pos_dict = {}
    amplicon = ''

    with open(f) as pos_file:
        for line in pos_file:
            if line.startswith('>'):
                amplicon = line.strip().lstrip('>')
                # can't have repeated headings, fine
                pos_dict[amplicon] = []
            else:
                line_split = line.strip().split(',')
                line_split[0] = int(line_split[0])
                line_split[1] = int(line_split[1])
                tup = (int(line_split[0]), int(line_split[1]), line_split[2])
                pos_dict[amplicon].append(line_split)

    return pos_dict

def add_fillbetween(ax, locs, color_dict=None):
    '''
    Adds highlighting to the region(s) defined in locs
    Locs should be a list of (start,end) tuples defining regions of interest
    We should have a dict somewhere where we can store motif names such that each only gets inserted into the legend once
    '''
    ymin, ymax = ax.get_ylim()

    # for (start, end, label, *color) in locs:
    for loc in locs:
        start, end, label, *color = loc

        if color_dict and label in color_dict:
            c = color_dict[label]
        else:
            c = color[0] if color else 'r'

        ax.fill_between((start,end), (ymin,ymin), (ymax,ymax), color=c, alpha=0.2, label=label)

def cluster_single_reads(mat):
    '''
    Cluster single reads in mat and return newly ordered matrix
    Note, this is super clunky since we are calling clustermap! No reason
    We shoudl actually just figure out whatever linkage it uses and return that (can do later)
    Also can optionally add in cols to use for clustering? 
    I don't think that any distance function will be tricked by entirely identical columns though...
    '''
    Z = linkage(mat)
    # cluster = fcluster(Z, criterion='maxclust', t=50)
    reorder = dendrogram(Z, no_plot=True)['leaves']
    return mat.iloc[reorder]
    
def plot_single_reads(mat, ax, fill_whitespace=False, axis_label=None):
    '''
    Plot single reads (clustered however) from mat on ax
    Mat should have 0 = methylated = accessible and 1 = unmethylated = protected
    -1 for positions where there is no data
    Then we map the colors and use gray cmap so 1->0=black and 0->0.7=gray (and -1->1=white)
    So black corresponds to protection and gray to open (and white to nada)
    Can optionally add back in the whitespace (if we've already removed those columns)
    '''
    cols = mat.columns
    if fill_whitespace:
        cols = list(range(1,max(mat.columns)))
    new_mat = pd.DataFrame(1-np.finfo(float).eps, index=mat.index, columns=cols)
    for c in mat.columns:
        # new_mat[c] = mat[c].map({1:0,0:0.7,-1:1-np.finfo(float).eps})
        new_mat[c] = mat[c].map({1:0,0:0.7,-1:1-np.finfo(float).eps})
        new_mat[c+1] = new_mat[c]
    ax.imshow(new_mat, aspect='auto', cmap=plt.get_cmap('gray'), interpolation='none')
    if axis_label:
        ax.set_ylabel(axis_label)

def fix_missing_data(mat):
    '''
    Substitute -1 (= missing) in single molecule matrix with 0.5 for clustering stuff
    '''
    return mat.replace(-1,0.5)

def parse_fasta(reference):
    """
    Parses the reference fasta file, and reads it into a dictionary. For the header, will
    strip anything after the first space. For example:
    
    >chr1 AC:XXX gi:XXX LN:XXX
    
    will be saved to the dictionary as:
    chr1
    """

    fasta_dict = {}
    fasta = open(reference, "r")
    for line in fasta:

        if line.startswith(">"): # Header line
            header = line.strip().split(" ")[0][1:] # Remove '>' and anything after ' '
            fasta_dict[header] = []

        else: # Sequence line
            fasta_dict[header].append(line.strip())
    fasta.close()    

    # If it's a multiline fasta, join the individual lines to one sequence
    for header in fasta_dict:
        fasta_dict[header] = "".join(fasta_dict[header])

    return fasta_dict

def get_methyl_positions(mat):
    '''
    Extracts meaningful columns from a single molecule methylation matrix
    Uninformative positions are stored as -1, so take only the columns with non-negative means 
    This isn't exactly correct but will do the trick for now, I think.

    FIX ME: Should instead compare 0/1 count to -1 count?
    '''

    return mat.loc[:,mat.mean()>=0].columns.tolist()

def load_single_molecule_matrix(path, name=None, every_other=False):
    '''
    Reads a single-molecule matrix from file, extracts the meaningful columns (using the function above) and returns
    Optionally (default) drops the duplicated columns (idk why but Georgi's script has 2 columns per GpC, one per C I guess?)
    Also, normally indexes molecules by their row number, however we often want to concatenate multiple rows, gives the option to prepend a unique sample name
    '''
    mat = pd.read_table(path, index_col=0)
    mat.columns = [int(c) for c in mat.columns] # convert cols to ints (helpful later)
    methyl_positions = get_methyl_positions(mat)
    if every_other:
        # filter out the duplicates (so each GpC only has one column)
        methyl_positions = [p for idx,p in enumerate(methyl_positions) if idx % 2 == 0]
    if name:
        mat.index = ['{}-{}'.format(name,i) for i in mat.index]
    return mat[methyl_positions].copy()

def load_single_molecule_classification(path, name=None):
    '''
    Reads a single-molecule classification matrix from file and returns it processed like the matrix above 
    Optionally add metadata columns (num bound and such)
    '''
    df = pd.read_table(path, index_col=0)
    if name:
        df.index = ['{}-{}'.format(name,i) for i in df.index]
    return df.copy()

def annotate_single_molecule_classifications(df):
    pass

def closest_gc(pos, gc_lst):
    '''
    Helper function which returns the nearest GpC for a given feature (pos)
    '''
    distances = [abs(pos-x) for x in gc_lst]
    return gc_lst[distances.index(min(distances))]

def plot_single_read(mat, idx, ax, label=True, fillbetween=None):
    '''
    Plot a single molecule trace on the axes provided
    Input is a matrix (mat) and the index of the molecule to plot (idx)
    Can optionally provide fillbetween = list of position tuples to add coloring
    1 is protected and 0 is accessible
    '''
    mat.loc[idx].reset_index().plot('index', idx, ax=ax, style='.-', c='k', ms=10)

    if fillbetween:
        add_fillbetween(ax, fillbetween)

def filter_all_converted_reads(mat, threshold=1):
    '''
    Function to take in a single-molecule matrix and remove all of the reads that are "entirely converted" = "entirely protected"
    We think these are artifacts, likely not totally heterochromatin? but rather conversion artifacts. Not sure though
    Pass in a single-molecule matrix and a threshold for what fraction of total bases need to be "protected" to qualify
    '''
    rows_to_exclude = mat.mean(axis=1) > threshold
    return mat.loc[~rows_to_exclude].copy()

def decorate_single_read_plot(assignment, ax, tfbs_positions):
    '''
    Takes an assignment of a single molecule and adds state calls to the plots
    State must be in the form of a dictionary with 'nuc_present' = boolean, 'nuc_pos' = (nuc_start, nuc_end), and tfbs_state = [list of booleans]
    Needs some coercing from the data table we output (see below)
    '''
    for idx,b in enumerate(assignment['tfbs_state']):
        if b:
            ax.add_patch(Ellipse(((0.5*tfbs_positions[idx][0]+0.5*tfbs_positions[idx][1]+3),0.1),10,0.1,color='r'))
    if assignment['nuc1_present']:
        ax.add_patch(Ellipse((0.5*assignment['nuc1_pos'][0]+0.5*assignment['nuc1_pos'][1],0.1),1*assignment['nuc1_pos'][1]-1*assignment['nuc1_pos'][0],0.1,color='0.4'))
    if assignment['nuc2_present']:
        ax.add_patch(Ellipse((0.5*assignment['nuc2_pos'][0]+0.5*assignment['nuc2_pos'][1],0.1),1*assignment['nuc2_pos'][1]-1*assignment['nuc2_pos'][0],0.1,color='0.7'))

def decorate_single_read_plot2(assignment, ax, tfbs_positions):
    '''
    
    '''
    # print(assignment)
    attributes = assignment.index.tolist()
    # print(attributes)

    # print(tfbs_positions)

    for idx in range(60):
        if 'nuc{}_present'.format(idx) in attributes and assignment['nuc{}_present'.format(idx)] == True:
            assignment['nuc{}_start'.format(idx)]
            ax.add_patch(Ellipse((0.5*assignment['nuc{}_start'.format(idx)]+0.5*assignment['nuc{}_end'.format(idx)],0.1),1*assignment['nuc{}_end'.format(idx)]-1*assignment['nuc{}_start'.format(idx)],0.1,color='0.4'))

    for idx in range(10):
        if 'tfbs_{}'.format(idx) in attributes and assignment['tfbs_{}'.format(idx)] == True:
            # TFBS indexed starting at 1 for some ungodly reason, so in order to get the corresponding one out of the positions file have to shift down the index
            ax.add_patch(Ellipse(((0.5*tfbs_positions[idx-1][0]+0.5*tfbs_positions[idx-1][1]+3),0.1),10,0.1,color='r'))

    promoter_components = {'TBP': (165,170,'darkblue'), 'PIC': (115,144,'lawngreen'), 'Pol2_Pause': (82,112,'gold')}

    for pc in promoter_components.keys():
        if pc in attributes and assignment[pc] == True:
            ax.add_patch(Ellipse(((0.5*promoter_components[pc][0]+0.5*promoter_components[pc][1]),0.1),1*promoter_components[pc][1]-1*promoter_components[pc][0],0.1,color=promoter_components[pc][2]))

def coerce_pandas_row_to_state_assignment(classification_row, n_tfbs): 
    '''
    Takes one row from the dictionary output by the bayesian classifier and converts it into the format we need to plot on single molecule plot
    '''
    state = {'nuc_present': classification_row['nuc_present'], 
             'nuc_pos': (classification_row['nuc_start'], classification_row['nuc_end']), 
             'tfbs_state': [classification_row['tfbs_{}'.format(i+1)] for i in range(n_tfbs)]}

    return state

def coerce_classification_dict_to_state_assignment(classification_dict, n_tfbs): 
    '''
    Takes one element dict from the immediate output of the bayesian pipeline (e.g. not yet written to disk)
    and convert to format for sm plotting
    Ok... it ends up being the same... but might not necessarily be in the future?
    '''
    state = {'nuc1_present': classification_dict['nuc1_present'], 
             'nuc1_pos': (classification_dict['nuc1_start'], classification_dict['nuc1_end']), 
             'nuc2_present': classification_dict['nuc2_present'], 
             'nuc2_pos': (classification_dict['nuc2_start'], classification_dict['nuc2_end']), 
             'tfbs_state': [classification_dict['tfbs_{}'.format(i+1)] for i in range(n_tfbs)]}

    return state

def adjust_gcgs(mat, gcg_pos, fix_pos):
    '''
    some sites are gcg (or cgc on the other strand or whatever)
    these sites aren't informative if methylated
    for some that are really close to others, we will just use the neighboring info to impute the methylated ones
    '''
    molecules_to_fix = mat[gcg_pos] == 0
    mat.loc[molecules_to_fix, gcg_pos] = mat.loc[molecules_to_fix, fix_pos]
    
    return mat

