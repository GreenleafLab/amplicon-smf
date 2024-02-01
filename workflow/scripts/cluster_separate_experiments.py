# input two matrices, plot clustered matrices separately
# then do KMeans on them together (and optionally separately) and return the cluster proportions for each
# change to barplot not pie chart
# then plot shapes of the clusters

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

from common import cluster_single_reads, plot_single_reads, fix_missing_data, add_fillbetween, load_tfbs_positions


def cluster_together(mat1, mat2, plots, methyl_positions, name1, name2, positions_file, amplicon, n_clusters, cmap=plt.get_cmap("tab10")):
    # do some general formatting
    colors = cmap([_ for _ in range(n_clusters)])
    title = amplicon if amplicon else ''
    name1 = name1 if name1 else '0'
    name2 = name2 if name2 else '1'

    # step1: idependently do hierarchical clustering on the reads and plot them adjacent
    mat1c, mat2c = cluster_single_reads(mat1), cluster_single_reads(mat2) # have to do this out of order because of clf with clustermap CHANGE

    fig, axs = plt.subplots(2,1,figsize=(12,10))
    ax_arr = axs.ravel()

    plot_single_reads(mat1c, ax_arr[0], name1)
    plot_single_reads(mat2c, ax_arr[1], name2)

    fig.suptitle('{} (Hierarchical)'.format(title))
    plt.xlabel('Position on plasmid')
    plots.savefig()
    plt.clf()

    # step 1.5 cluster all the reads with KMeans together
    # potentially abstract this out into common.py so we can do this also on single experiments
    df1 = mat1.copy()
    df2 = mat2.copy()
    df1['sample'] = name1
    df2['sample'] = name2
    df12 = df1.append(df2)
    labels = df12['sample']
    df12 = df12.drop('sample', axis=1)
    kmeans = KMeans(n_clusters=n_clusters).fit(df12)
    df12['cluster'] = kmeans.labels_
    df12['sample'] = labels

    # step2: jointly do KMeans on all the reads and plot the reads jointly organized by cluster
    # NOTE: would be nice here to use sns clustermap so we can add colorbars to the sides indicating which cluster they belong to...
    fig, axs = plt.subplots(2,1, figsize=(12,10))
    ax_arr = axs.ravel()

    plot_single_reads(df12.loc[df12['sample']==name1].sort_values('cluster').drop(['sample','cluster'], axis=1), ax_arr[0], name1)
    plot_single_reads(df12.loc[df12['sample']==name2].sort_values('cluster').drop(['sample','cluster'], axis=1), ax_arr[1], name2)
    
    fig.suptitle('{} (Joint KMeans)'.format(title))
    plt.xlabel('Position on plasmid')
    plots.savefig()
    plt.clf()

    # step 3: plot the cluster proportions used by each sample
    fig, ax = plt.subplots(figsize=(12,8))

    # make sure colors are the same between
    pd.crosstab(df12['sample'], df12['cluster']).apply(lambda r: r/r.sum()*100, axis=1).plot.bar(ax=ax, stacked=True)#, colors=colors)

    ax.set_ylabel('Cluster proportions')
    ax.set_xlabel('Sample')

    fig.suptitle('{} (Cluster Usage)'.format(title))
    plots.savefig()
    plt.clf()

    # step 4: plot the cluster shapes on a graph together
    df12 = fix_missing_data(df12)
    
    fig, ax = plt.subplots(figsize=(10,6))
    for i in range(n_clusters):
        ax.plot(methyl_positions, df12.loc[df12['cluster']==i,methyl_positions].mean()*100, color=colors[i], label=i)

    pos_dict = None
    if positions_file:
        pos_dict = load_tfbs_positions(positions_file)
    if pos_dict and amplicon:
        # have to write something here to make sure that the amplicon name georgi style and my style are the same...
        positions = pos_dict[amplicon]
        add_fillbetween(ax, positions)

    ax.set_ylabel('%SMF')
    ax.set_ylim((0,103))

    plt.legend()
    fig.suptitle('{} (Cluster shapes)'.format(title))
    plt.xlabel('Position on plasmid')
    plots.savefig()
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--input1", dest="input1", type=str, help="Path to single molecule matrix1")
    parser.add_argument("--input2", dest="input2", type=str, help="Path to single molecule matrix2")
    # parser.add_argument("--output", dest="output", type=str, help="Prefix for writing the merged output tables to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--positions", dest="positions_file", type=str, default=None, help="Optional path to file with TFBS sites to color")
    parser.add_argument("--input1_name", dest="input1_name", type=str, default=None, help="Name of input1 sample (e.g. -dox)")
    parser.add_argument("--input2_name", dest="input2_name", type=str, default=None, help="Name of input1 sample (e.g. +dox)")
    parser.add_argument("--amplicon_name", dest="amplicon_name", type=str, default=None, help="Name of amplicon aligned to (for plots and also to grab TFBS positions)")
    parser.add_argument("--n_clusters", dest="n_clusters", type=int, default=6, help="Number of clusters to use for kmeans")
    parser.add_argument("--reads_to_use", dest="reads_to_use", type=int, default=0, help="Number of reads to sample for clustering")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        # read files
        mat1 = pd.read_table(args.input1, index_col=0)
        mat1.columns = [int(x) for x in mat1.columns]
        mat2 = pd.read_table(args.input2, index_col=0)
        mat2.columns = [int(x) for x in mat2.columns]

        methyl_positions = mat1.loc[:,mat1.mean()>0].columns.tolist()
        # methyl_positions2 = mat2.loc[:,mat2.mean()>0].columns.tolist()
        # assert methyl_positions = methyl_positions2

        if args.reads_to_use:
            print(args.reads_to_use)
            print(mat1.head())
            mat1 = mat1.sample(n=args.reads_to_use).copy()
            mat1.index = list(range(args.reads_to_use))
            print(mat1.head())
            mat2 = mat2.sample(n=args.reads_to_use).copy()#.reset_index(inplace=False)
            mat2.index = list(range(args.reads_to_use))


        cluster_together(mat1, mat2, plots, methyl_positions, args.input1_name, args.input2_name, args.positions_file, args.amplicon_name, args.n_clusters)


