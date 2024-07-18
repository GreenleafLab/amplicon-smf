##################################
#                                #
# Last modified 2018/04/24       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import string
import pysam
import numpy as np
import random
import os
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sets import Set
import pandas as pd

# FLAG field meaning
# 0x0001 1 the read is paired in sequencing, no matter whether it is mapped in a pair
# 0x0002 2 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
# 0x0004 4 the query sequence itself is unmapped
# 0x0008 8 the mate is unmapped 1
# 0x0010 16 strand of the query (0 for forward; 1 for reverse strand)
# 0x0020 32 strand of the mate 1
# 0x0040 64 the read is the first read in a pair 1,2
# 0x0080 128 the read is the second read in a pair 1,2
# 0x0100 256 the alignment is not primary (a read having split hits may have multiple primary alignment records)
# 0x0200 512 the read fails platform/vendor quality checks
# 0x0400 1024 the read is either a PCR duplicate or an optical duplicate

def getReverseComplement(preliminarysequence):
    
    DNA = {'A':'T','T':'A','G':'C','C':'G','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    sequence=''
    for j in range(len(preliminarysequence)):
        sequence=sequence+DNA[preliminarysequence[len(preliminarysequence)-j-1]]
    return sequence

def FLAG(FLAG):

    Numbers = [0,1,2,4,8,16,32,64,128,256,512,1024]

    FLAGList=[]

    MaxNumberList=[]
    for i in Numbers:
        if i <= FLAG:
            MaxNumberList.append(i)

    Residual=FLAG
    maxPos = len(MaxNumberList)-1

    while Residual > 0:
        if MaxNumberList[maxPos] <= Residual:
            Residual = Residual - MaxNumberList[maxPos]
            FLAGList.append(MaxNumberList[maxPos])
            maxPos-=1
        else:
            maxPos-=1
  
    return FLAGList

def run():

    if len(sys.argv) < 9:
        print 'usage: python %s BAM genome.fa CG|GC|both peak_list chrFieldID leftFieldID rightFieldID strandFieldID outfile_prefix [-subset N] [-label fieldID] [-noEndogenousMethylation] [-minCov fraction] [-unstranded] [-heatmap path_to_heatmap.py x_pixel_size y_pixel_size colorscheme width(inches,dpi)] [-cluster]' % sys.argv[0]
        print '\tThe BAM file should come from Bismarck, but is expected to be sorted and deduped'
        print '\tUse the [-subset] option if you want only N of the fragments; the script will pick the N fragments best covering each region, and will discard regions with fewer than N covering fragments'
        print '\tUse the [-label] option if you want regions to be labeled with something other than their coordinates'
        print '\tBy default the script will filter out CpGpC sites if only [GpC] was given as the methylation type; use the [-noEndogenousMethylation] option to included those too'
        print '\tThe [-heatmap] option will generate png heatmaps instead of text file matrices'
        print '\tThe [-minCov] option will remove all fragments that cover the region at less than the specified fraction'
        print '\tNote: the script assumes no indels in the BAM file!!!'
        print '\tNote: the script assumes PBAT libraries by default!!!'
        print '\tNote on output: '
        print '\t\tscore of 0 means methylated GC/CG, i.e. no protection'
        print '\t\tscore of 1 means unmethylated GC/CG, i.e. footprinting protection'
        print '\t\tscore of -1 means no GC or CG'
        sys.exit(1)

    BAM = sys.argv[1]
    fasta = sys.argv[2]
    GCtype = sys.argv[3]
    peaks = sys.argv[4]
    chrFieldID = int(sys.argv[5])
    leftFieldID = int(sys.argv[6])
    rightFieldID = int(sys.argv[7])
    strandFieldID = int(sys.argv[8])
    outprefix = sys.argv[9]
    outdonefile = sys.argv[10]
#     amplicon = sys.argv[10]

    GenomeDict={}
    sequence=''
    inputdatafile = open(fasta)
    for line in inputdatafile:
        if line[0]=='>':
            if sequence != '':
                GenomeDict[chr] = ''.join(sequence).upper()
            chr = line.strip().split('>')[1]
            sequence=[]
            Keep=False
            continue
        else:
            sequence.append(line.strip())
    GenomeDict[chr] = ''.join(sequence).upper()

    GCDict={}
    for chrom in GenomeDict.keys():
        seq = GenomeDict[chrom]
        gc_positions = [i for i,b in enumerate(seq) if b=='G' and seq[min(i+1,len(seq)-1)]=='C']
        GCDict[chrom] = gc_positions
    # print GCDict

    GCGDict={}
    for chrom in GenomeDict.keys():
        seq = GenomeDict[chrom]
        gcg_positions = [i for i,b in enumerate(seq) if b=='G' and seq[min(i+1,len(seq)-1)]=='C' and seq[min(i+2,len(seq)-1)]=='G']
        GCGDict[chrom] = gcg_positions
    # print GCGDict['opJS2_1xTetO']

    print 'finished inputting genomic sequence'

    doNoEndMeth = False
    if '-noEndogenousMethylation' in sys.argv:
        doNoEndMeth = True

    doSubset = False
    if '-subset' in sys.argv:
        doSubset = True
        Nsub = int(sys.argv[sys.argv.index('-subset')+1])
        print 'will only output the most complete', Nsub, 'fragments'

    doMincov = False
    if '-minCov' in sys.argv:
        doMincov = True
        MFC = float(sys.argv[sys.argv.index('-minCov')+1])
        print 'will only output fragments with at least', MFC, 'fractional coverage of the region'

    doLabel = False
    if '-label' in sys.argv:
        doLabel = True
        labelFieldID = int(sys.argv[sys.argv.index('-label')+1])

    doHM = False
    if '-heatmap' in sys.argv:
        doHM = True
        print 'will output heatmap'
        HMpy = sys.argv[sys.argv.index('-heatmap')+1]
        HMxp = int(sys.argv[sys.argv.index('-heatmap')+2])
        HMyp = int(sys.argv[sys.argv.index('-heatmap')+3])
        HMcs = sys.argv[sys.argv.index('-heatmap')+4]
        HMincdpi = sys.argv[sys.argv.index('-heatmap')+5]

    doNS = False
    if '-unstranded' in sys.argv:
        doNS = True
        print 'will not treat regions as stranded'

    doCluster = False
    if '-cluster' in sys.argv:    
        doCluster = True
        print 'will attempt to cluster the reads'

    samfile = pysam.Samfile(BAM, "rb" )

    with open(outdonefile, 'w') as out:
        out.write('amplicon\ttotal_reads\tobserved_states\treads_per_state\n')

        linelist = open(peaks)
        for line in linelist:
            print(line)
            if line.startswith('#'):
                continue
            # if amplicon not in line:
            #     continue
            # print(amplicon)
            linefields = line.strip().split('\t')
            chr = linefields[chrFieldID]
            start = max(0,int(linefields[leftFieldID]))
            end = int(linefields[rightFieldID])
            if doNS:
                strand = '+'
            else:
                strand = linefields[strandFieldID]
            ReadDict = {}
            Matrix = []
            if doLabel:
                label = linefields[labelFieldID]
            else:
                if strand == '+':
                    label = chr + '_' + str(start) + '-' + str(end) + '_for'
                if strand == '-':
                    label = chr + '_' + str(start) + '-' + str(end) + '_rev'
            for alignedread in samfile.fetch(chr, start, end):

                fields=str(alignedread).split('\t')
                ID = fields[0]
                if ReadDict.has_key(ID):
                    pass
                else:
                    ReadDict[ID] = []
                FLAGfields = FLAG(int(fields[1]))
                pos = alignedread.pos - 1
    #            readseq = alignedread.seq
                readseq_temp = alignedread.seq
                readseq = ''
                rpos = 0
                for (m,bp) in alignedread.cigar:
    # soft-clipped bases:
                    if m == 4:
                        rpos += bp
    # matches:
                    if m == 0:
                        readseq += readseq_temp[rpos:rpos+bp]
                        rpos += bp
    # insertions:
    # note: not handled properly, as the junction remaining after excising the insertion might be a CG or GC
    # but there is no good way to deal with it
                    if m == 1:
                        rpos += bp
    # deletions:
                    if m == 2:
                        for D in range(bp):
                            readseq += 'N'
                ReadDict[ID].append((FLAGfields,pos,readseq))
    #            if alignedread.is_reverse:
    #                s = '-'
    #            else:
    #                s = '+'
    #            MD = alignedread.opt('MD')
    #            if alignedread.is_read2:
    #                ReadDict[ID][2] = 
    #            if alignedread.is_read1:
    #                ReadDict[ID][1] = 
            for ID in ReadDict.keys():
                # print('---')
                # print(ID)
                scores = []
                for i in range(start,end):
                    scores.append(-1)
                covered = []
                for (FLAGfields,pos,readseq) in ReadDict[ID]:
                    # print(readseq)
    #                genomeseq = GenomeDict[chr][pos:pos + len(readseq)]
    #                if (128 in FLAGfields and 16 in FLAGfields) or (64 in FLAGfields and 32 in FLAGfields):
    #                if (128 in FLAGfields and 32 in FLAGfields) or (64 in FLAGfields and 16 in FLAGfields):
                    for i in range(pos,pos + len(readseq)):
                        if i < start or i >= end:
                            continue
                        covered.append(i)
                        if GCtype == 'both':
                            if GenomeDict[chr][i:i+2] == 'CG':
                                if readseq[i-pos-1:i-pos+1] == 'CG':
                                    scores[i-start] = 0
                                    scores[min(i+1-start,end-start-1)] = 0
                                else:
                                    scores[i-start] = 1
                                    scores[min(i+1-start,end-start-1)] = 1
                            if GenomeDict[chr][i:i+2] == 'GC':
                                if readseq[i-pos-1:i-pos+1] == 'GC':
                                    scores[i-start] = 0
                                    scores[min(i+1-start,end-start-1)] = 0
                                else:
                                    scores[i-start] = 1
                                    scores[min(i+1-start,end-start-1)] = 1
                        if GCtype == 'GC':
                            if GenomeDict[chr][i:i+2] == 'GC':
                                if doNoEndMeth:
                                    pass
                                else:
                                    # if GenomeDict[chr][i-1:i+2] == 'CGC':
                                    if GenomeDict[chr][i:i+3] == 'GCG':
                                        continue
                                if readseq[i-pos-1:i-pos+1] == 'GC':
                                    scores[i-start] = 0
                                    scores[min(i+1-start,end-start-1)] = 0
                                else:
                                    scores[i-start] = 1
                                    scores[min(i+1-start,end-start-1)] = 1
                        if GCtype == 'CG':
                            if GenomeDict[chr][i:i+2] == 'CG':
                                if readseq[i-pos-1:i-pos+1] == 'CG':
                                    scores[i-start] = 0
                                    scores[min(i+1-start,end-start-1)] = 0
                                else:
                                    scores[i-start] = 1
                                    scores[min(i+1-start,end-start-1)] = 1

                    # problematic_positions = 0
                    # for gc_pos in GCDict[chr]:
                    #     if gc_pos not in GCGDict[chr]:
                    #         if gc_pos > pos and gc_pos < pos + len(readseq):
                    #             if scores[gc_pos-start:gc_pos-start+2] == [-1, -1]:
                    #                 problematic_positions += 1
                    #                 print gc_pos
                    #                 print GenomeDict[chr][gc_pos-2:gc_pos+2]
                    #                 print readseq[gc_pos-pos-2:gc_pos-pos+2]
                    # print(problematic_positions)

                covered = Set(covered)
                if doMincov:
                    if len(covered)/(end-start-0.0) < MFC:
                        continue
                Matrix.append((len(covered),scores,ID)) # append the readID as well so include in the unclustered matrix output
            if len(Matrix) < 1: 
                print('No reads for: ' + label)
                continue
            if doSubset:
                if len(Matrix) < Nsub:
                    NewMatrix = Matrix
                else:
                    Matrix.sort()
                    Matrix.reverse()
                    if Matrix[Nsub-1][0] < Matrix[0][0]:
                        NewMatrix = Matrix[0:Nsub]
                    else:
                        for i in range(len(Matrix)):
                            if Matrix[i][0] < Matrix[0][0]:
                                TempMatrix = Matrix[0:i]
                                break
                            TempMatrix = Matrix[0:i]
                        print len(TempMatrix), Nsub
                        if Nsub > len(TempMatrix):
                            continue
                        NewMatrix = random.sample(TempMatrix,Nsub)
                print 'M, NM', len(Matrix), len(Matrix[0][1]), len(NewMatrix), len(NewMatrix[0][1]) 
            else:
                NewMatrix = Matrix
            print len(NewMatrix), 'reads retained for', label
            if len(NewMatrix) <= 1:
                continue

            # here is where we edit Georgi script
            # we will write out the entirety of Matrix to a file, in the order in which we saw the reads
            # then we will optionally try to do the clustering thing below, on NewMatrix

            # reuse the output code for saving the full matrix
            outfile = open(outprefix + '.' + label + '.full_unclustered.matrix', 'w')
            outline = '#' + chr
            for i in range(start,end):
                outline = outline + '\t' + str(i)
            outfile.write(outline+'\n')
            for idx, (L,scores,ID) in enumerate(Matrix):
                if strand == '-':
                    scores.reverse()
                outline = str(ID)                   # write the read ID as the first column value instead of the redudant "idx"
                for s in scores:
                    outline = outline + '\t' + str(s)
                outfile.write(outline + '\n')
            outfile.close()

            print 'wrote full matrix to file: ' + outprefix + '.' + label + '.full_unclustered.matrix'

            # if we want to cluster, do the clustering on NewMatrix
            if doCluster:
                X = []
                for (L,scores,ID) in NewMatrix:
                    X.append(scores)
                Z = linkage(X, method='ward', metric='euclidean', optimal_ordering=True)
                clusters = fcluster(Z, 0, criterion='distance')
                CDict = {}
                for i in range(len(clusters)):
                    C = clusters[i]
                    if CDict.has_key(C):
                        pass
                    else:
                        CDict[C] = []
                    CDict[C].append(i)
                Cs = CDict.keys()
                Cs.sort()

                outfile = open(outprefix + '.' + label + '.clustered.matrix', 'w')
                outline = '#' + chr
                for i in range(start,end):
                    outline = outline + '\t' + str(i)
                outfile.write(outline+'\n')
                for C in Cs:
                    for k in CDict[C]:
                        scores = X[k]
                        if strand == '-':
                            scores.reverse()
                        outline = str(C)
                        for s in scores:
                            outline = outline + '\t' + str(s)
                        outfile.write(outline + '\n')
                outfile.close()
                print 'wrote clustered matrix to file: ' + outprefix + '.' + label + '.clustered.matrix'

        #        sys.exit(1)
                if doHM:
                    cmd = 'python ' + HMpy + ' ' + outprefix + '.' + label + '.clustered.matrix' + ' ' + str(HMxp) + ' ' + str(HMyp) + ' -0.25 1 ' + HMcs + ' ' + HMincdpi + ' ' + outprefix + '.' + label + '.matrix.png'
                    contents = os.system(cmd)
                    # cmd = 'rm ' + outprefix + '.' + label + '.matrix'
                    # contents = os.system(cmd)

            # now we need to output some QC metrics on the reads, e.g. what fraction of the total reads have unique methylation states
            df = pd.DataFrame([x[1] for x in Matrix])
            gc_cols = df.loc[:,df.mean()>=0].columns.tolist()
            df = df[gc_cols]

            total_reads = len(df)
            unique_methyl_states = len(df.drop_duplicates())
            reads_per_state = (total_reads+0.0) / unique_methyl_states

            out.write('{}\t{}\t{}\t{:.2f}\n'.format(label,total_reads,unique_methyl_states,reads_per_state))

    # with open(outdonefile, 'w') as out:
    #     out.write('Done clustering')
            
run()


