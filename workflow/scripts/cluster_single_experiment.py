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


def cluster(mat, plots, methyl_positions, name, positions_file, amplicon, n_clusters, cmap=plt.get_cmap("tab10")):
    # do some general formatting
    colors = cmap([_ for _ in range(n_clusters)])
    title = amplicon if amplicon else ''
    name = name if name else ''

    # step1: idependently do hierarchical clustering on the reads and plot them adjacent
    matc = cluster_single_reads(mat) # have to do this out of order because of clf with clustermap CHANGE

    fig, ax = plt.subplots(1,1,figsize=(12,5))

    plot_single_reads(matc, ax, name)

    plt.title('{} (Hierarchical)'.format(title))
    plt.xlabel('Position on plasmid')
    plots.savefig()
    plt.clf()

    # step 1.5 cluster all the reads with KMeans together
    # potentially abstract this out into common.py so we can do this also on single experiments
    df = mat.copy()
    kmeans = KMeans(n_clusters=n_clusters).fit(df)
    df['cluster'] = kmeans.labels_

    # step2: jointly do KMeans on all the reads and plot the reads jointly organized by cluster
    # NOTE: would be nice here to use sns clustermap so we can add colorbars to the sides indicating which cluster they belong to...
    fig, ax = plt.subplots(1,1, figsize=(12,5))

    plot_single_reads(df.sort_values('cluster').drop('cluster', axis=1), ax, name)
    
    plt.title('{} (KMeans)'.format(title))
    plt.xlabel('Position on plasmid')
    plots.savefig()
    plt.clf()

    # step 3: plot the cluster proportions used by each sample
    fig, ax = plt.subplots(figsize=(12,8))

    # make sure colors are the same between
    pd.DataFrame((df['cluster'].value_counts(normalize=True)*100)).T.plot.bar(ax=ax, stacked=True)
    # pd.crosstab(df12['sample'], df12['cluster']).apply(lambda r: r/r.sum()*100, axis=1).plot.bar(ax=ax, stacked=True)#, colors=colors)

    ax.set_ylabel('Cluster proportions')
    ax.set_xlabel('Sample')

    plt.title('{} (Cluster Usage)'.format(title))
    plots.savefig()
    plt.clf()

    # step 4: plot the cluster shapes on a graph together
    df = fix_missing_data(df)
    
    fig, ax = plt.subplots(figsize=(10,6))
    for i in range(n_clusters):
        ax.plot(methyl_positions, df.loc[df['cluster']==i,methyl_positions].mean()*100, color=colors[i], label=i)

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
    plt.title('{} (Cluster shapes)'.format(title))
    plt.xlabel('Position on plasmid')
    plots.savefig()
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plots bulk methylation signal across amplicons')
    parser.add_argument("--input", dest="input", type=str, help="Path to single molecule matrix")
    # parser.add_argument("--output", dest="output", type=str, help="Prefix for writing the merged output tables to")
    parser.add_argument("--plot", dest="plot_file", type=str, help="Path to output plots file")
    parser.add_argument("--positions", dest="positions_file", type=str, default=None, help="Optional path to file with TFBS sites to color")
    parser.add_argument("--input_name", dest="input_name", type=str, default=None, help="Name of input sample (e.g. CTCF)")
    parser.add_argument("--amplicon_name", dest="amplicon_name", type=str, default=None, help="Name of amplicon aligned to (for plots and also to grab TFBS positions)")
    parser.add_argument("--n_clusters", dest="n_clusters", type=int, default=6, help="Number of clusters to use for kmeans")

    args = parser.parse_args()

    with PdfPages(args.plot_file) as plots:
        # read files
        mat = pd.read_table(args.input, index_col=0)
        mat.columns = [int(x) for x in mat.columns]

        methyl_positions = mat.loc[:,mat.mean()>0].columns.tolist()

        cluster(mat, plots, methyl_positions, args.input_name, args.positions_file, args.amplicon_name, args.n_clusters)


