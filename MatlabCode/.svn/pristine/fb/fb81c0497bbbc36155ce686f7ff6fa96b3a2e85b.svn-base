# Code to help me find candidate promoter sequences that we can put into AAV
# Contents:
#
#
# Section 1
# Getting a count of gene hits for all Magno and Parvo, integrating across all of our RNA-seq
# datasets. Then loading in MACS2 output for two ATAC-seq datasets. This code will need to be
# changed if we ever want to pool across > 1 ATAC-seq dataset.

import os
import numpy as np
from matplotlib import pyplot
import HTSeq
import copy

RNAseqpathroot = '/Volumes/2TBdrive/RNA'
RNAseqsubpaths = ['RNA_2016_01', 'RNA_2016_05']
RNAseqfilenames = ['A09177_M_q_genecounts', 'A09177_P_q_genecounts',
                   'A20416_P2_q_genecounts', 'A21015_M_q_genecounts',
                   'A20416_M_q_genecounts','A20416_P1_q_genecounts']
Magnofiles = {'A09177_M_q_genecounts','A21015_M_q_genecounts','A20416_M_q_genecounts'}
Parvofiles = {'A09177_P_q_genecounts','A20416_P1_q_genecounts','A20416_P2_q_genecounts'}

# Making a dictionary of gene names as keys and lists as values where the first element in each
# list is the count from the first RNAseqfile and the second element is the count from the second
RNAdictlist = []
for filename in RNAseqfilenames:
    file_in = None
    for subpath in RNAseqsubpaths:
        wholefilename = os.path.join(RNAseqpathroot, subpath, filename)
        if os.path.isfile(wholefilename):
            file_in = open(wholefilename)
    if file_in is None:
        raise IOError("Couldn't find file: {0}".format(wholefilename))
    tmpRNAdict = {}
    while True:
        line = file_in.readline()
        linelist = line.split()
        if (line == '' or linelist[0][0] == '_'):
            break
        assert linelist[0] not in tmpRNAdict, "two genes with the same name detected"
        tmpRNAdict[linelist[0]] = int(linelist[1])
    file_in.close()
    RNAdictlist.append(tmpRNAdict) # a list of dictionaries; one for each RNA-seq file.

RNAdict = {}
for key in RNAdictlist[0].keys():
    tmpcounts = np.zeros(shape=(1,2),dtype = int)
    for i in range(len(RNAdictlist)):
        if RNAseqfilenames[i] in Magnofiles: # somehow I have to remember that 0 = magno and 1 = parvo
           tmpcounts[0,0] += RNAdictlist[i][key]
        else:
           tmpcounts[0,1] += RNAdictlist[i][key]
        del RNAdictlist[i][key]
    RNAdict[key] = tmpcounts[:] # <--- Somehow this [:] makes this a deep copy which I need

for i in RNAdictlist: assert not i, "Extra genes detected in one or more RNA dictionaries."
RNAarray = np.zeros(shape=(len(RNAdict),2), dtype = int)
gene_names = RNAdict.keys()
for i in range(len(gene_names)): RNAarray[i,:] = RNAdict[gene_names[i]]

# At this point the RNA-seq data is represented in two ways:
# 1) RNAdict: Keys are genes and values are 1x2 ndarrays of counts
# 2) gene_names: the keys are RNAdict + RNAarray: nx2 array of counts
#         in the same order as gene_names
# IMPORTANT: "gene_names" is going to be the gold standard for now on.
# that will help me keep things in order.

RNATHRESHOLD_ABS = .0001
RNATHRESHOLD_REL = .1
normalizedRNAarray = np.transpose(np.vstack((RNAarray[:,0]/float(sum(RNAarray[:,0])),RNAarray[:,1]/float(sum(RNAarray[:,1])))))
pyplot.plot(normalizedRNAarray[:,0],normalizedRNAarray[:,1],'.')
pyplot.axis('equal')
L = ((normalizedRNAarray[:,0] > RNATHRESHOLD_ABS) & (normalizedRNAarray[:,1] > RNATHRESHOLD_ABS))
logdiffs = np.log(normalizedRNAarray[:,0]) - np.log(normalizedRNAarray[:,1])
L = L & (abs(logdiffs) > RNATHRESHOLD_REL)
pyplot.plot(normalizedRNAarray[L,0],normalizedRNAarray[L,1],'r.')
[gene_names[i] for i,j in enumerate(L) if j] # seems to be working
LRNA = L

# Now pulling in the ATAC-seq data
# I guess the idea here is to start with a dictionary of peak tuples (start, stop, fold increase, p-value) that are
# in the vicinity of each TSS. How do I process these data into a form that will identify big, differential peaks?
# I guess I cam go through the dict coming up with some kind of score per TSS.

# First finding the TSSs
# SLOW.I should save this in a ile somewhere.
gtffile = HTSeq.GFF_Reader("/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf")
tssdict = {}
for feature in gtffile:
    if feature.type == "exon" and feature.attr["exon_number"] == "1" :
        p = feature.iv.start_d_as_pos  # start_d_as_pos = The position of the TSS irrespective of strand
        #if (p.start > halfwinwidth):
        #    window = HTSeq.GenomicInterval(p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, feature.iv.strand)
        #    tsspos[window] += (feature.name, p.pos)
        tssdict[feature.name] = [p.chrom, p.pos]

ATACseqpathroot = '/Volumes/2TBdrive/ATAC'
ATACseqsubpaths = ['ATAC_2016_05','ATAC_2016_09']
ATACseqfilenames = ['Z12131_M_q_peaks_peaks.narrowPeak','Z12131_P_q_peaks_peaks.narrowPeak']

# Collecting all the peaks from each ATACseqfile, creating a GenomicArrayofSets (like for the TSSs)
# and saving them in a list.
# I think the order of columns in the macs2 output file are ... fold change, -log10(p), -log10(q), step level (for browser)
ATAClist = [] # list of GenomicArrayofSets (one for M and one for P)
for filename in ATACseqfilenames:
    file_in = None
    for subpath in ATACseqsubpaths:
        wholefilename = os.path.join(ATACseqpathroot, subpath, filename)
        if os.path.isfile(wholefilename):
            file_in = open(wholefilename)
    if file_in is None:
        raise IOError("Couldn't find file: {0}".format(wholefilename))
    tmpATACdict = {} # keys are chroms. values are 4-tuples: start, stop, q, fold increase
    # NOTE: online description of MACS2 output file format is wrong. Also the space between
    # "2TB" and "drive" was a problem here.
    while True:
        line = file_in.readline()
        linelist = line.split()
        if (line == '' or linelist[0][0] == '_'):
            break
        if (linelist[0] not in tmpATACdict):
            tmpATACdict[linelist[0]] = [(int(linelist[1]), int(linelist[2]), float(linelist[8]), float(linelist[6])),]
        else:
            tmpATACdict[linelist[0]].append((int(linelist[1]), int(linelist[2]), float(linelist[8]), float(linelist[6])))
    file_in.close()
    # Now moving the contents of the dict to a GenomicArrayOfSets
    peakpos = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    for chrom in tmpATACdict:
        for v in tmpATACdict[chrom]:
            window = HTSeq.GenomicInterval(chrom, v[0], v[1]) # start and stop
            peakpos[window] += (v[2], v[3]) # -log10(q-value) and fold increase
    ATAClist.append(copy.deepcopy(peakpos))

# For testing the above code
iv = HTSeq.GenomicInterval('chr8',1,70000000)
print(list(ATAClist[0][iv].steps()))
print(list(ATAClist[1][iv].steps()))

#############################
# Going to go through each TSS and finding peaks in the neighborhood.
# Stats is the sum of peak fold increases in the vicinity of each TSS.
# (Which might not even really be the fold increase if I have the columns wrong)
#############################
halfwinwidth = 3000
stats = np.zeros(shape=RNAarray.shape)
for k in gene_names:
    iv = HTSeq.GenomicInterval(tssdict[k][0], max(0, tssdict[k][1] - halfwinwidth), tssdict[k][1] + halfwinwidth)
    try:
        for ATACdictidx in range(len(ATAClist)):  # Looping over M and P
            for v in ATAClist[ATACdictidx][iv].steps():  # looping over chromosomes
                for i in v[1]:
                    stats[gene_names.index(k), ATACdictidx] += i[1]
    except KeyError:
        continue

# Trying to find signifcant peaks in one that doesn't overlap with a signifcant peak in the other.
halfwinwidth = 1000
stats = np.zeros(shape=RNAarray.shape)
for k in gene_names: # Looping over TSSs
    iv = HTSeq.GenomicInterval(tssdict[k][0], max(0, tssdict[k][1] - halfwinwidth),
                                   tssdict[k][1] + halfwinwidth) # window around each TSS
    s0 = list()  # initializing
    s1 = list()  # initializing
    for ATACidx in range(len(ATAClist)):  # Looping over M and P
        try:
            for v,s in ATAClist[ATACidx][iv].steps():  # finding peaks near TSS
                if s: # If s is "None" there's no peak over the interval v
                    if ATACidx == 0:
                        s0.append(v)
                    else:
                        s1.append(v)
        except KeyError:
            continue
    #if s0 == s1:
     #    print('peaks are the same')
     #    break

    if s0 or s1:
        print('got here')
        # Counting the number of significant peaks that are within halfwinwidth of TSS
        # and do not have an overlapping peak in the other sample.
        n_unique_peaks = [[1]*len(s0), [1]*len(s1)]
        for i in s0:
            for j in s1:
                if (i.overlaps(j)):
                    n_unique_peaks[0][s0.index(i)] = 0
                    n_unique_peaks[1][s1.index(j)] = 0
        stats[gene_names.index(k), :] = np.asarray([sum(n_unique_peaks[0]),sum(n_unique_peaks[1])])

# Testing to see whether this is analysis is bringing our attention
# to peaks near TSS

pyplot.plot(stats[:,0],stats[:,1],'.')
pyplot.plot(stats[LRNA,0],stats[LRNA,1],'r.')

L = LRNA & (stats[:,0] == 0) & (stats[:,1]> 1.5)
L = LRNA & (stats[:,0] > 1.5) & (stats[:,1] == 0)
#L = LRNA & (stats[:,0] > 10) & (stats[:,1]> 10)
print([gene_names[i] for i,j in enumerate(L) if L[i]])
