##################################
#                                #
# Last modified 2021/01/11       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import os
import string
import math
import numpy as np
import matplotlib, copy
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.patches import Circle, RegularPolygon
from matplotlib.collections import PatchCollection

def run():

    if len(sys.argv) < 8:
        print 'usage: python %s datafile x_pixel_size y_pixel_size min|min max|max colorscheme width(inches,dpi) outfile [-every Nth] [-top N] [-average rows columns] [-saveps] [-RowLabels] [-ColLabels]' % sys.argv[0]
        print '\tNote: the input file can be .bz2 or gz'
        print '\tNote: the script assumes that the left most column contains entry IDs'
        print '\tNote: enter width as a comma-separated tuple of the inches and dpi; the height will be rescaled accordingly'
        print '\tNote: http://matplotlib.org/examples/color/colormaps_reference.html'
        print '\tNote: the [-top N] option assumes the entries have been already sorted'
        sys.exit(1)
    
    input = sys.argv[1]
    xps = float(sys.argv[2])
    yps = float(sys.argv[3])
    if sys.argv[4] != 'min':
       minS = float(sys.argv[4])
    else:
       minS = 'min'
    if sys.argv[5] == 'max':
       maxS = 'max'
    else:
       maxS = float(sys.argv[5])
    cscheme = sys.argv[6]
    (inches,DP) = sys.argv[7].split(',')
    inches = float(inches)
    DP = int(DP)
    outfilename = sys.argv[8]

    doPostScript = False
    if '-saveps' in sys.argv:
        doPostScript = True

    doSkip = False
    if '-every' in sys.argv:
        doSkip = True
        SkipInterval = int(sys.argv[sys.argv.index('-every') + 1])
        print 'will only show every', SkipInterval, 'entry'

    doTop = False
    if '-top' in sys.argv:
        doTop = True
        TopN = int(sys.argv[sys.argv.index('-top') + 1])
        print 'will only show the top', TopN, 'entries'

    doXlabels = False
    if '-RowLabels' in sys.argv:
        doXlabels = True
        print 'will label rows'
        XLs = []

    doYlabels = False
    if '-ColLabels' in sys.argv:
        doYlabels = True
        print 'will label columns'
        YLs = []

    DataMatrix = []

    if minS == 'min' or maxS == 'max':
        if input.endswith('.bz2'):
            cmd = 'bzip2 -cd ' + input
        elif input.endswith('.gz'):
            cmd = 'zcat ' + input
        else:
            cmd = 'cat ' + input
        p = os.popen(cmd, "r")
        line = 'line'
        i = 0
        while line != '':
            line = p.readline()
            if line == '':
                break
            fields = line.strip().split('\t')
            if line.startswith('#'):
                continue
            i += 1
            if i % 10000 == 0:
                print i
            if doSkip:
                if i % SkipInterval != 0:
                    continue
            if i == 1:
                NewMin = float(fields[1])
                NewMax = float(fields[1])
            for j in range(1,len(fields)):
                v = float(fields[j])
                NewMin = min(NewMin,v)
                NewMax = max(NewMax,v)
        if minS == 'min':
            minS = NewMin
        if maxS == 'max':
            maxS = NewMax
        print minS, maxS

    if input.endswith('.bz2'):
        cmd = 'bzip2 -cd ' + input
    elif input.endswith('.gz'):
        cmd = 'zcat ' + input
    else:
        cmd = 'cat ' + input
    p = os.popen(cmd, "r")
    line = 'line'
    i = 0
    while line != '':
        line = p.readline()
        if line == '':
            break
        fields = line.strip().split('\t')
        if line.startswith('#'):
            if doYlabels:
                for j in range(1,len(fields)):
                    YLs.append(fields[j])
            continue
        i+=1
        if i % 10000 == 0:
            print i
        if doSkip:
            if i % SkipInterval != 0:
                continue
        row = []
        if doXlabels:
            XLs.append(fields[0])
        for j in range(1,len(fields)):
            if fields[j] == 'nan':
                v = 0
            else:
                v = float(fields[j])
            if v < minS:
                v = minS
            if v > maxS:
                v = maxS
            v = (v - minS)/(maxS-minS)
            row.append(v)
        DataMatrix.append(row)

    if doTop:
        DataMatrix = DataMatrix[0:TopN]

    DataMatrix.reverse()
    if doXlabels:
        XLs.reverse()


    NRows = len(DataMatrix)
    NColumns = len(DataMatrix[0])

    Height = NRows*yps
    Width = NColumns*xps

    print NRows, yps, Height
    print NColumns, xps, Width
    print inches, Height/Width, inches*(Height/Width)

    if '-average' in sys.argv:
        AX = int(sys.argv[sys.argv.index('-average') + 1])
        AY = int(sys.argv[sys.argv.index('-average') + 2])
        NewDataMatrix = []
        for i in range(0,len(DataMatrix),AX):
            x_array = []
            for j in range(0,len(DataMatrix[0]),AY):
                i_array = []
                for k1 in range(i,min(i + AX,len(DataMatrix))):
                     for k2 in range(j,min(j + AY,len(DataMatrix[0]))):
                         i_array.append(DataMatrix[k1][k2])
                x_array.append(np.mean(i_array))
            NewDataMatrix.append(x_array)

        DataMatrix = NewDataMatrix

        NRows = len(DataMatrix)
        NColumns = len(DataMatrix[0])

        Height = NRows*yps
        Width = NColumns*xps

        print NRows, yps, Height
        print NColumns, xps, Width
        print inches, Height/Width, inches*(Height/Width)

    print len(DataMatrix)
    print len(DataMatrix[0])

    if doXlabels or doYlabels:
        rect = 0.09,0.05,0.9,0.9
    else:
        rect = 0,0,1,1
    # fig = figure(figsize=(80, 20),dpi=100)
    fig = figure(figsize=(inches, inches*(Height/Width)),dpi=DP)
    ax = fig.add_subplot(1,1,1,aspect='equal')
    ax = fig.add_axes(rect)

#    ax.pcolor(DataMatrix, cmap=plt.cm.Blues, alpha=0.8)
#    ax.imshow(np.array(DataMatrix), vmin=0, vmax=1, cmap=cscheme, aspect = 'auto')
    ax.pcolor(DataMatrix, vmin=0, vmax=1, cmap=cscheme)

    ax.set_xticks([])
    ax.set_yticks([])
    if doXlabels:
        ax.set_yticks(np.arange(+.5, len(XLs), 1))
#        ax.set_ylabel(XLs, fontsize=xps)
        ax.set_yticklabels(XLs)
    if doYlabels:
        ax.set_xticks(np.arange(+.5, len(YLs), 1))
#        ax.set_xlabel(YLs, fontsize=yps)
        ax.set_xticklabels(YLs)

#    ax.set_yticks(np.arange(len(XLs)))
#    ax.set_yticklabels(XLs)

#    ax.grid(which='major', axis = 'Y', color='k', linestyle='-', linewidth=2)

    if outfilename.endswith('.png'):
        savefig(outfilename)
    else:
        savefig(outfilename + '.png')

    if doPostScript:
        if outfilename.endswith('.png'):
            savefig(outfilename[0:-4] + '.eps', format='eps')
        else:
            savefig(outfilename + '.eps', format='eps')
   
run()
