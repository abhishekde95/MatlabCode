# This is going to be a script for analyzing ATAC-seq reads near NUMTs.
# Starting by seeing whether we can make TSS-aligned histograms segregrated on the basis of RNA expression.

import sys
import HTSeq
from matplotlib import pyplot
import numpy as np
import itertools
import os
import dill
import pysam
import MapReduceHackery as MRH

# Making a TSS plot using MapReduce
bamfilename1 = MRH.findfile('A13063_V1_40min_sorted.bam');
bamfilename1 = MRH.findfile('Z13145_Str_q_sorted.bam');

gapfilename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/gap.txt' 
blacklistfilename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/rhemac2blacklist.bed'
halfwinwidth = 3000

# Loads the "tsslist"
try:
    workspace = dill.load(open('/Users/horwitzlab/Desktop/MatlabCode/Greg/PythonStuff/workspace.p','r'))
except:
    print("Cannot find 'workspace'. Need to run this commented out code, below.")
    sys.exit()
    # Getting all the TSSs in a list (tsslist).
    # Commenting this out. Only need this to recalculate tsslist.
    # Takes about a couple of minutes to run.
    # making tsslist a list of user-defined objects doesn't work because the parallel processing stuff
    # does not work with "non-pickleable" objects
    # tsslist = []
    # tsspos = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    # gtffile = HTSeq.GFF_Reader("/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf")
    # for feature in gtffile:
    #    if feature.type == "exon" and feature.attr["exon_number"] == "1":
    #       p = feature.iv.start_d_as_pos  # start_d_as_pos = The position of the TSS irrespective of strand
    #       if (p.pos - halfwinwidth > 0):
    #           tsslist.append([p.chrom,  p.strand, p.pos, feature.name]) # chrom, strand, tsspos, name
    # tsslist.append(tss(p.chrom,  p.strand, p.pos, feature.name)) # chrom, strand, tsspos, name
    # dill.dump(tsslist, open('workspace.p','w'))
tsslist = workspace
tsslist = [tsslist[i] for i in range(len(tsslist)) if tsslist[i][0][0:3] == 'chr' and tsslist[i][0] != 'chrMT'] # Getting rid of TSSs on the mitochondrial genome or the random contigs

gaplesstsslist = MRH.EliminateGapsFromTssList(tsslist, halfwinwidth, gapfilename)
print('Eliminating {} TSSs of {} ({}%) because of gaps.'.format(str(len(tsslist)-len(gaplesstsslist)), str(len(tsslist)),
                                                                str(100*float(len(tsslist)-len(gaplesstsslist))/float(len(tsslist)))))
    
profile = MRH.CalculateTSSHist(bamfilename1, goodchromsonlytsslist, halfwinwidth, blackListFileName = blacklistfilename)
pyplot.figure()
pyplot.plot(range(-halfwinwidth,halfwinwidth),profile)
pyplot.title(bamfilename1[bamfilename1.rfind(os.path.sep)+1:])
pyplot.show()

#############################################################
#
# As above, but looping through genes based on expression level
#
#############################################################

RNAseqpath = '/Volumes/2TBdrive/RNA/RNA_2016_05'
RNAseqfilename = 'A20416_M_q_genecounts'
try:
    file = open(RNAseqpath+os.path.sep+RNAseqfilename)
except IOError as inst:
    print(RNAseqpath+os.path.sep+RNAseqfilename+' does not seem to exist')
else:
    RNAdict = {}
    while True:
        line = file.readline()
        linelist = line.split()
        if line == '':
            break
        elif linelist[0] not in RNAdict and linelist[0][0] != '_':
            RNAdict[linelist[0]] = int(linelist[1])
        elif linelist[0] in RNAdict:
            print("two genes with the same name detected")
            break
file.close()
maxreads = max(RNAdict.values())
genehitprctiles = [25, 50, 60, 70, 80, 90]
cutoffvalues = np.percentile(RNAdict.values(),genehitprctiles)
cutoffvalues = np.append(cutoffvalues, maxreads)

# Getting the ranks of the RNA expression levels (using raw read counts)
RNAranks = {}
for a, b in RNAdict.iteritems():
    for j in range(len(cutoffvalues)):
       if a not in RNAranks and b <= cutoffvalues[j] :
           RNAranks[a] = j # ranks of gene expression
           break

gaplesstsslist = MRH.EliminateGapsFromTssList(tsslist, halfwinwidth)
profiles = np.zeros(shape=(2*halfwinwidth,len(cutoffvalues)))
numberofgenes = [0]*len(cutoffvalues)
for expressionrankcounter in range(len(cutoffvalues)): # This would go faster if we weren't making n passes through the GTF file
    print('expressionrankcounter is {0}'.format(expressionrankcounter))
    subtsslist = [i for i in gaplesstsslist if RNAranks[i[3]] == expressionrankcounter]
    profiles[:,expressionrankcounter] = MRH.CalculateTSSHist(bamfilename1, subtsslist, halfwinwidth, blackListFileName = blacklistfilename)
    numberofgenes[expressionrankcounter] = len(subtsslist)

for i in range(len(numberofgenes)):
   pyplot.plot(np.arange(-halfwinwidth, halfwinwidth),profiles[:,i]/float(numberofgenes[i]))

legendstr = ['< '+str(genehitprctiles[0])]
if len(genehitprctiles) > 2:
   legendstr += [str(genehitprctiles[i])+'-'+str(genehitprctiles[i+1])+'%' for i in range(len(genehitprctiles)-1)]
legendstr += ['> '+str(genehitprctiles[-1])]
pyplot.legend(legendstr,fontsize = 10)
pyplot.ylabel('reads per TSS')
pyplot.xlabel('Distance from TSS (bp)')
pyplot.title(bamfilename1[bamfilename1.rfind(os.path.sep)+1:])


###########################################################
#
# Aligning ATAC-seq reads on the begining and ends of NUMTs
#
############################################################
NUMTs_filename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/NUMTs.txt'
#NUMTs_filename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/NUMTs_all.bed'
file_in = open(NUMTs_filename,'rb')
NUMTs_list = []
while True:
    line = file_in.readline()
    linelist = line.split()
    if (line == '' or linelist[0][0] == '_'):
        break
    if NUMTs_filename == '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/NUMTs.txt':
       NUMTcenter = int(np.mean([int(linelist[2]), int(linelist[3])])) # if using NUMTs.txt
       NUMTs_list.append([linelist[0], linelist[1], NUMTcenter])
    else:
        NUMTcenter = int(np.mean([int(linelist[1]), int(linelist[2])])) # if using numts-uniq.bed
        if (len(linelist) == 4):
            NUMTs_list.append([linelist[0], linelist[4], NUMTcenter])
        else:
            NUMTs_list.append([linelist[0], NUMTcenter])

file_in.close()

profile = MRH.CalculateTSSHist(bamfilename1, NUMTs_list, halfwinwidth)
pyplot.plot(np.arange(-halfwinwidth, halfwinwidth),profile)
pyplot.title('ATAC-seq reads aligned to NUMTs')
pyplot.xlabel('Distance from NUMTs center')

# Now doing the counting
# need to recalculate NUMTs_list because we need start and stop separately, not their mean
file_in = open(NUMTs_filename,'rb')
NUMTs_list = []
while True:
    line = file_in.readline()
    linelist = line.split()
    if (line == '' or linelist[0][0] == '_'):
        break
    NUMTs_list.append([linelist[0], linelist[1], int(linelist[2]), int(linelist[3])])
    #NUMTs_list.append([linelist[0], linelist[4], int(linelist[1]), int(linelist[2])])
file_in.close()
NumberOfReadsOverlappingNUMTs = MRH.CountBAMOverlapReads(bamfilename1,NUMTs_list)

# --------------------------------------------------
# How much of the genome is occupied by these NUMTs?
totalNUMTarea = 0
for i in NUMTs_list:
    totalNUMTarea = totalNUMTarea + (i[3]-i[2])
# How big is the genome?
bam = pysam.AlignmentFile(bamfilename1, "rb")
GenomeLength = sum(bam.lengths)
chrom_sizes = zip(bam.references, bam.lengths)
bam.close()
# Fraction of genome covered by NUMTs
totalNUMTarea/float(GenomeLength)

# --------------------------------------------------
# Making a version of NUMTs_list that spans the entire genome
whole_genome_list = [[i[0], '.', 0,i[1]] for i in chrom_sizes]
NumberOfReadsOverall = MRH.CountBAMOverlapReads(bamfilename1,whole_genome_list) # takes a while but using the same algorithm
NumberOfReadsOverlappingNUMTs/float(NumberOfReadsOverall)

#%%
# Counting up how many reads map onto which chromosome
import os
import numpy as np
from matplotlib import pyplot
from __future__ import division


REVERSEFILEORDER = True # Put recent files in the front. Important since we get the chromosome names from the first bamfile processed
samtoolsstr = "/Users/horwitzlab/Desktop/SEQ_ANALYSIS/samtools-1.3\ 2/samtools"
ATACmaindir = '/Volumes/2TBdrive/ATAC/ATAC_2017_02'
#bamfilenames = ['Z13145_Str_q_sorted.bam', 'Z13145_V1t_q_sorted.bam', 'Z13145_V1b_q_sorted.bam', 'Z12131_M_sorted.bam','Z12131_P_sorted.bam','Z12131_V1t_sorted.bam', 'Z12131_V1b_sorted.bam']# requires index
bamfilenames = ['A13062_V1_nuc_sorted.bam','A13063_V1_10min_sorted.bam','A13063_V1_40min_sorted.bam']
tmpsubdirs = !{"ls "+ATACmaindir}

# Getting the full paths for the bam files and reordering the contents of bamfilenames to be consistent with "bamfiles"
bamfiles = []
tmpfilenames = []
for i in range(len(tmpsubdirs)):
    if os.path.isdir(ATACmaindir+os.sep+tmpsubdirs[i]):
        for filename in bamfilenames:
            fullpath = ATACmaindir+os.sep+tmpsubdirs[i]+os.sep+filename
            if os.path.isfile(fullpath):
                bamfiles.append(fullpath)
                tmpfilenames.append(filename)
bamfilenames = tmpfilenames
if REVERSEFILEORDER:
    bamfiles.reverse()
    bamfilenames.reverse()
    
# Makes a quick pass through the first bamfile to get the chromosome names.
# Also initializing the readcounts array.
# If subsequent bamfiles have yet other chromosomes, ignoring them.
samtools_out = !{samtoolsstr + " idxstats " + bamfiles[0]}
chrnames = []
chrlengths = []
for line in samtools_out:
    tmp = line.split()
    if (tmp[0][0:3] == 'chr'):
        chrnames.append(tmp[0])
        chrlengths.append(tmp[1])
chrlengths = np.array(chrlengths)# changing data type     
readcounts = np.empty((len(chrnames),len(bamfiles)))
readcounts[:] = np.NAN    
# Now getting read numbers. Each bamfile is a column, each chromosome is a row.
for bamfileidx in range(len(bamfiles)):
    bamfile = bamfiles[bamfileidx]
    samtools_out = !{samtoolsstr + " idxstats " + bamfile}
    for line in samtools_out:
        tmp = line.split()
        idx = [i for i,j in enumerate(chrnames) if tmp[0] == j]
        if (idx):
            readcounts[idx,bamfileidx] = tmp[2]
pyplot.figure()
pyplot.plot(readcounts)
pyplot.legend(bamfilenames,fontsize = 7,loc = 0)
pyplot.ylabel('# reads')
abbreviatedchrnames = [i[3:] for i in chrnames]
pyplot.xticks(range(len(abbreviatedchrnames)),abbreviatedchrnames)
pyplot.xlabel('chromosome')

MTidx = [i for i in range(len(chrnames)) if chrnames[i] == 'chrMT'] # getting the index for the chrMT
for i in range(len(tmpfilenames)):
    print(tmpfilenames[i]+': {0} which is {1}% of the total').format(readcounts[MTidx,i], 100*readcounts[MTidx,i]/sum(readcounts[:,i]))

# Getting proportions of reads per chromosome
pyplot.figure()
for i in range(readcounts.shape[1]):
    pyplot.plot(readcounts[:,i]/chrlengths.astype(float))
pyplot.yscale('log')
pyplot.ylabel('readcounts/chrom length')
abbreviatedchrnames = [i[3:] for i in chrnames]
pyplot.xticks(range(len(abbreviatedchrnames)),abbreviatedchrnames)
pyplot.xlabel('chromosome')
pyplot.legend(bamfilenames,fontsize = 7,loc = 0)
#%%
# Histogram of read lengths
import MapReduceHackery as MRH
import deeptools.mapReduce
import pysam
import numpy as np

blacklistfilename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/rhemac2blacklist.bed'
bamfilename = MRH.findfile('Z12131_V1b.bam');
readcounts = MRH.CalculateReadLengthHist(bamfilename, maxreadlength = 3000, blackListFileName = blacklistfilename)
pyplot.figure()
pyplot.plot(readcounts)
pyplot.title(bamfilename[bamfilename.rfind(os.path.sep)+1:])
