# Playing around with the DeepTools API. Working through https://media.readthedocs.org/pdf/deeptools/documentation/deeptools.pdf
# For what it's worth, I think I'm using deeptools version 2.2.4
#
#  8/23/16 GDLH

import deeptools.countReadsPerBin as countReadsPerBin
import os
import matplotlib.pyplot as plt
import numpy as np

# import matplotlib.axis.Axis as axis

ATACseqpath = '/Volumes/2TBdrive/ATAC/ATAC_2017_02/Horwitz'
ATACseq_bam_files = ('A13062_V1_nuc_sorted.bam','A13063_V1_10min_sorted.bam','A13063_V1_10min_sorted.bam')
#bam_files = ('Z12125_M_sorted.bam',)
# CountReadsPerBin needs an index, so the reads in the BAM file must be sorted by position, not name
bam_file_list = list()
for filename in ATACseq_bam_files:
    bam_file_list.append(ATACseqpath+os.path.sep+filename)

# FINDING COVERAGE OVER A REGION
cr = countReadsPerBin.CountReadsPerBin(bam_file_list, binLength=10, stepSize=10)
# cr is a "CountReadsPerBin" object
start_end_pos = (120000, 221000)
readtuple = cr.count_reads_in_region('chr2',start_end_pos[0],start_end_pos[1])
# docs say the output of cr.count_reads_in_region shold be an ndarray, but it's a tuple. First element is
# an ndarray like I was hoping for, but I have no idea what the second element (an empty string) is for.
# An older version, 1.6, returned a tuple. Am I using 1.6?!
readarray = readtuple[0]
plt.plot(np.arange(start_end_pos[0],start_end_pos[1],cr.binLength),readarray)

# SAMPLING THE GENOME
totalSites = 10000
cr = countReadsPerBin.CountReadsPerBin(bam_file_list, binLength=1, numberOfSamples=totalSites,
                                        numberOfProcessors=8) # The skipZeros argument isn't working
sequencing_depth = cr.run()
print sequencing_depth.mean(axis=0) # number of reads overlapping the average base in the genome

# Trying the plotting in the documentation
fig, axs = plt.subplots(1, 2, figsize=(15,5))
# plot coverage
for col in sequencing_depth.T: # iterating over the number of bamfiles
    axs[0].plot(np.bincount(col.astype(int)).astype(float)/totalSites) # bascially doing "hist"
    csum = np.bincount(col.astype(int))[::-1].cumsum()
    axs[1].plot(csum.astype(float)[::-1] / csum.max())
axs[0].set_xlabel('coverage')
axs[0].set_ylabel('fraction of bases sampled')
# plot cumulative coverage

axs[1].set_xlabel('coverage')
axs[1].set_ylabel('fraction of bases sampled >= coverage')

# ---------------------------------------------------
# Comparing read depth between MiSeq and HighSeq runs
seqpath = '/Volumes/2TBdrive/ATAC'
seq_bam_files = ('ATAC_2016_05/Z12125_P2_sorted.bam','ATAC_2016_09/Z12131_p_sorted.bam')
bam_file_list = list()
for filename in seq_bam_files:
    bam_file_list.append(seqpath+os.path.sep+filename)
cr = countReadsPerBin.CountReadsPerBin(bam_file_list, binLength=1, numberOfSamples=totalSites,
                                        numberOfProcessors=8)
sequencing_depth = cr.run()
print sequencing_depth.mean(axis=0)

# Trying to make a log-spaced histogram just for kicks.
bins = np.logspace(-1,2.0, num = 50)
plt.hist(sequencing_depth, bins)
plt.gca().set_xscale('log')

#####################################################
## Just trying to import each module from deeptools
# which I now realize is a "package", not a module
# "SES" = "signal extraction scaling" from Diaz et al. (2012)
import deeptools.SES_scaleFactor
test = deeptools.SES_scaleFactor.Tester()
bin_length = 50000
num_samples = 400
with Timer():
   _dict = deeptools.SES_scaleFactor.estimateScaleFactor([bamfilename1, bamfilename2],
                                                      bin_length, num_samples,  1, numberOfProcessors = 8)
print(_dict['size_factors'])
print(_dict['size_factors_based_on_mean'])

# %%
######################################################
#
from __future__ import division
import os
import operator
import numpy as np
from matplotlib import pylab as plt
import deeptools.bamCompare as bamCompare   # little help
import deeptools.bamHandler as bamHandler   # no help
import deeptools.bigwigCompare as bigwigCompare   # no help
import deeptools.computeGCBias as computeGCBias  # lots of help
import deeptools.computeMatrix as computeMatrix  # no help
import deeptools.plotHeatmap as plotHeatmap
import deeptools.multiBamSummary as multiBamSummary
import deeptools.plotCorrelation as plotCorrelation
import deeptools.plotPCA as plotPCA
#
##########################################################
# Trying out bamCoverage -> computeMatrix -> plotHeatmap/plotProfile
##########################################################
import deeptools.bamCoverage as bamCoverage
BLACKOUTGAPS = False # Better be "False" unles you want to get stuck in a loop
#bamfilenamestem = 'Z12131_V1t_q_sorted'
bamfilenamestem = 'A13062_V1_nuc_sorted'
gapfilename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/gap.txt'
args = ['-b','/Volumes/2TBdrive/ATAC/ATAC_2017_02/Horwitz/'+bamfilenamestem+'.bam']
args = args + ['--outFileFormat', 'bigwig', '-bs', '1', '-p', '7','-o','coverage.bw','--samFlagInclude','66'] # 66 = first read, read in proper pair

if BLACKOUTGAPS: # doesn't work! System hangs
    args += ['-bl',gapfilename]
os.chdir('/Users/horwitzlab/Desktop')
bamCoverage.main(args) # creates a file called "coverage.bw" 

# Now calling computeMatrix
# Works with GTF file but not BED file
import deeptools.computeMatrix as computeMatrix
gtffilename = "/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf"
outfilename = 'computeMatrixOut_'+bamfilenamestem # some binary format
args = ['reference-point','-S','coverage.bw','-R',gtffilename,'-a','3000','-b','3000','-out', outfilename,'--referencePoint','TSS','-p','7']
if BLACKOUTGAPS:
    args += ['-bl',gapfilename]
computeMatrix.main(args)

# Now calling plotHeatmap
import deeptools.plotHeatmap as plotHeatmap
args = ['-m',outfilename,'-out','plotHeatMapOut','--verbose']
plotHeatmap.main(args)

# Now calling plotProfile
import deeptools.plotProfile as plotProfile
args = ['-m',outfilename,'-out','plotProfileOut','--verbose']
plotProfile.main(args)

#%%
####################
# Trying bamCompare
# Just making coverage vectors from bam files and seeing where in the genome they differ.
# This worked but wasn't super useful
####################
#  The "args" input to main should be a list
ATACseqpath = '/Volumes/2TBdrive/ATAC/ATAC_2017_02/Horwitz'
ATACseq_bam_files = ('A13063_V1_10min_sorted.bam','A13063_V1_40min_sorted.bam','A13063_V1_nuc_sorted')

arg  = ['-b1', ATACseqpath+os.path.sep+ATACseq_bam_files[0], '-b2', ATACseqpath+os.path.sep+ATACseq_bam_files[1],
        '-bs','100000','-o' ,'bamCompareOut','-of','bedgraph']
greg = bamCompare.main(arg) # The output is "None" is everything worked. The output file is a BedGraph

# Making two passes through the Bedgraph file - once to get the distribition of log2 differences and another to identify the regions
# of the genome where these differences occur.
GWdifferences = dict()
try:
    file = open('bamCompareOut')
except IOError as inst:
    print('The file bamCompareOut does not seem to exist')
else:
    while True: # This takes a while because we're going through the whole genome
        line = file.readline()
        if line == '':
            break
        linelist = line.split()
        if float(linelist[3]) in GWdifferences:  # keeping track of log2 read lengths?
            GWdifferences[float(linelist[3])] += 1
        else:
            GWdifferences[float(linelist[3])] = 1
    file.close()

sorted_GWdifferences = sorted(GWdifferences.items(), key=operator.itemgetter(0))
plt.figure()
for i in sorted_GWdifferences:
    plt.plot(i[0],i[1],'.')
plt.show()
threshold = 1
# Now making a second pass through the file and finding regions that differ by more than threshold
interesting_genome_regions = [];
try:
    file = open('bamCompareOut')
except IOError as inst:
    print('The file bamCompareOut does not seem to exist')
else:
    counter = 0
    while True: # This takes a while because we're going through the whole genome
        line = file.readline()
        if line == '':
            break
        counter += 1
        linelist = line.split()
        if abs(float(linelist[3])) > threshold:
            interesting_genome_regions.append(linelist)
        if (counter % 10000 == 0 ):
            print(counter)
    file.close()
interesting_genome_regions.sort()
for i in interesting_genome_regions:
    print(i[0]+':'+i[1]+'-'+i[2]+'         '+i[3])
#%%
##########################
# Trying multiBamSummary
# This makes the scatterplot of ATAC-seq coverage between samples
# This confirmed that there is a bimodal distribution of coverage
# and that the two modes are consistent between comparisons
##########################
ATACseqpathroot = ('/Volumes/2TBdrive/ATAC')
ATACseq_bam_files = ('A13063_V1_10min_sorted.bam','A13063_V1_40min_sorted.bam','A13062_V1_nuc_sorted.bam')
ATACseq_bam_files_withpath = []
for (dirpath, dirnames, filenames) in os.walk(ATACseqpathroot):
    for a in ATACseq_bam_files:
         if a in filenames:
             ATACseq_bam_files_withpath.extend([dirpath+os.path.sep+a])
    
NUMTS_file = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/rhemac2blacklist.bed'
#NUMTS_file = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/numts-uniq.bed'
BLACKOUTNUMTS = False
arg  = ['bins', '-b']
arg+=ATACseq_bam_files_withpath
arg+= ['-out', 'MBSout.npz','--ignoreDuplicates', '--numberOfProcessors','6','--outRawCounts','readCounts.tab','--samFlagInclude','66']
if BLACKOUTNUMTS:
    arg += ['-bl', NUMTS_file]
multiBamSummary.main(arg)

arg = ['--corData','MBSout.npz','-o','pltcorout','-c','spearman','--whatToPlot','heatmap']
plotCorrelation.main(arg) # saves a file
arg = ['--corData','MBSout.npz','-o','PCAout']
plotPCA.main(arg) # saves a file

sampleidxs = [0, 2] # Which samples to compare
if np.max(sampleidxs) > len(ATACseq_bam_files)-1:
    print('Sample index out of range')
else:
    greg = np.load('MBSout.npz')
    labels = greg['labels']
    matrix = greg['matrix']
    plt.figure()
    plt.plot(np.log10(matrix[:,sampleidxs[0]]),np.log10(matrix[:,sampleidxs[1]]),'r.')
    plt.plot([0,5], [0,5],'r-')
    plt.xlabel(ATACseq_bam_files[sampleidxs[0]])
    plt.ylabel(ATACseq_bam_files[sampleidxs[1]])

# Tracking down the positions of putative NUMTs
L = np.logical_and((np.log10(matrix[:,sampleidxs[0]]) < np.log10(matrix[:,sampleidxs[1]])-.2), (np.log10(matrix[:,sampleidxs[0]]) > 2)) # categories based on M vs P comparison
plt.plot(np.log10(matrix[L,sampleidxs[0]]),np.log10(matrix[L,sampleidxs[1]]),'b.')

# Looking for bimodality of the ratio of coverage
# It seems like it's the really big ATAC-seq peaks
# that differ the most between whatever these two cell groups are
thresholds = np.logspace(1, 6, 16)
plt.figure()
counter = 0
for i in range(len(thresholds)):
    ax = plt.subplot(np.ceil(np.sqrt(len(thresholds))), np.ceil(np.sqrt(len(thresholds))), counter+1)
    counter += 1
    if (i == 0):
        L_lower = np.logical_and(matrix[:,0] > 0, matrix[:,1] > 0)
    else:
        L_lower = np.logical_and(matrix[:,0] > thresholds[i-1], matrix[:,1] > thresholds[i-1])
    L_upper = np.logical_and(matrix[:,0] < thresholds[i], matrix[:,1] < thresholds[i])
    L = np.logical_and(L_lower, L_upper)
    plt.hist(np.log2(matrix[L,0])-np.log2(matrix[L,1]), 100)
    plt.xlim(-1, 1)
    if i < len(thresholds) :
        plt.xticks([])
    plt.title(str(np.round(thresholds[i])))

# Showing that the collection of open genomic regions is consistent
# across comparisons (M vs P, V1t vs V1b, etc)
L = matrix[:,1] > 5/7*matrix[:,0] # categories based on M vs P comparison
axislabels = ['M','P','V1t','V1b']
plt.figure()
counter = 1
for i in range(len(axislabels)):
    for j in range(i+1,len(axislabels)):
        ax = plt.subplot(2,3,counter)
        plt.plot(matrix[L,i],matrix[L,j],'b.')
        plt.plot(matrix[~L,i],matrix[~L,j],'r.')
        counter += 1
        ax.set_xlabel(axislabels[i])
        ax.set_ylabel(axislabels[j])
        plt.xticks([0, 250000])
        plt.yticks([0, 250000])
#%%
# -----------------------------
# Making a map of the genome with marks indicating where the big differential peaks are
# -----------------------------
chromsizefile = open('/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/RheMacChromSizes.txt')
chrsizedict = {}
while True:
    line = chromsizefile.readline()
    linelist = line.split()
    if line == '':
        break
    if linelist[0] == 'chrUr':
        continue
    chrsizedict[linelist[0]] = int(linelist[1])
chromsizefile.close()
sizevect=np.zeros(len(chrsizedict))
# in sizevect, elements 0-19 are chrs 1 through 20 and element 20 is chrX
for a in range(1,21):
    sizevect[a-1] = chrsizedict['chr'+str(a)]
sizevect[-1] = chrsizedict['chrX']
plt.figure()
for i in range(0,len(sizevect)):
    plt.plot([0, sizevect[i]],[i+1, i+1],'k-')
plt.ylim(0, 23)



#%%
# -----------------------------
# Looking at the distribution of insert sizes per chromosome. 
# Are the mitochondrial reads the shortest? (Can we get rid of them)
# using bead selection?
# -----------------------------
import MapReduceHackery as MRH
import numpy as np
from matplotlib import pyplot

chromsizefile = open('/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/RheMacChromSizes.txt')
chrsizelist = []
while True:
    line = chromsizefile.readline()
    linelist = line.split()
    if line == '':
        break
    if linelist[0] == 'chrUr':
        continue
    chrsizelist.append((linelist[0],int(linelist[1])))
chromsizefile.close()

bamfilename = MRH.findfile('A13063_V1_40min_sorted.bam')
data = np.array([])
for i in range(len(chrsizelist)):
    a = MRH.CalculateReadLengthHist (bamfilename, chrom_sizes = [chrsizelist[i]] , maxreadlength = 2000)
    data = np.vstack([data, a]) if data.size else a # Clever python code for concatenation of nparray!

pyplot.figure()
for i in range(len(chrsizelist)):
   if chrsizelist[i][0] == 'chrMT'
pyplot.plot(np.transpose(data))
pyplot.yscale('log')

#%%

# Finding regions of the chromosome where M and P both have lots of reads
# Note, the lines in readCounts.tab are not in the same order as the rows
# in "matrix". That's why I'm not using matrix, but overwriting it.
# readCounts.tab is generated by multiBamSummary.
coveragesummary = open('readCounts.tab')
line = coveragesummary.readline() # Getting rid of the header line
nreadsarray = np.zeros(shape=(matrix.shape[0],4)) # <--- this '4' needs to be changed
chrposarray = np.zeros(dtype='int',shape=(matrix.shape[0],2))
chrnames = [None for x in range(len(matrix))]
counter = 0
while True:
    line = coveragesummary.readline()
    linelist = line.split()
    if line == '':
        break
    chrposarray[counter,:] = linelist[1:3]
    nreadsarray[counter,:] = linelist[3:7]
    chrnames[counter] = linelist[0]
    counter += 1
coveragesummary.close()
L = np.logical_and(nreadsarray[:,1] > 5/7*nreadsarray[:,0], nreadsarray[:,1] > 10000)
print("P is dominant for these (nuclear open chromatin?)")
idxs = [i for i, x in enumerate(L) if x]
for i in idxs:
    print(str(chrnames[i])+':'+str(chrposarray[i,0])+'-'+str(chrposarray[i,1]))

L = np.logical_and(nreadsarray[:,1] < 5/7*nreadsarray[:,0], nreadsarray[:,1] > 10000)
print("M is dominant for these (mitochondrial open chromatin?)")
idxs = [i for i, x in enumerate(L) if x]
for i in idxs:
    print(str(chrnames[i])+':'+str(chrposarray[i,0])+'-'+str(chrposarray[i,1]))

L = np.logical_and(nreadsarray[:,0] > 10000, nreadsarray[:,1] > 10000)
for i, j in enumerate(L):
    if j:
        if chrnames[i] == 'chrX':
            ypos = 21
        elif chrnames[i] == 'chrUr':
            continue
        else:
            ypos = int(chrnames[i][3:])
        xpos = np.mean(chrposarray[i,0:1])
        h = plt.plot(xpos, ypos,'*')
        if nreadsarray[i, 1] > 5 / 7 * nreadsarray[i, 0]: # P is dominant
            h[0]._markeredgecolor = 'blue'
            h[0]._markerfacecolor = 'blue'
        else:
            h[0]._markeredgecolor = 'red'
            h[0]._markerfacecolor = 'red'
plt.title('Open chromatin by chromosome and cluster')
plt.xlabel('genomic postion (bp)')
plt.ylabel('chromosome')

########################################################################
# Looking at one of those coincident ATAC-seq/RNA-seq peaks and
# asking whether the reads are on the same strand or on different strands.
# Maybe the ATAC-seq peaks are due to mitochondrial DNA contamination
# and the RNA-seq peaks are due to mitochondrial RNA (or DNA?) contamination
# In any case, it'll be good to know that I can do this analysis.
# ########################################################################
import deeptools.bamCoverage as bamCoverage
whichchrom = 'chr2'
start = 174574711
end = 174574849

whichchrom = 'chr14'
start = 48503020
end = 48503526

whichchrom = 'chr16'
start = 55497127
end = 55507784

bamfilenamestems = ['Z12131_M_sorted','Z12131_P_sorted','A20416_M_sorted','A20416_P2_sorted']
ATACdir = '/Volumes/2TBdrive/ATAC/ATAC_2016_09'
RNAdir = '/Volumes/2TBdrive/RNA/RNA_2016_05'
data = np.zeros(shape=(end-start,2*len(bamfilenamestems)))
bedgraphfilename = 'tmpbedgraph'

for bfs in bamfilenamestems:
    if os.path.isfile(ATACdir+os.sep+bfs+'.bam'):
        args = ['-b',ATACdir+os.sep+bfs+'.bam']
    elif os.path.isfile(RNAdir+os.sep+bfs+'.bam'):
        args = ['-b',RNAdir+os.sep+bfs+'.bam']
    else:
        print("Cannot find file {} anywhere I looked.".format(bfs+'.bam'))
        continue
    filenameidx = bamfilenamestems.index(bfs)
    for direction in ['forward','reverse']:
        diridx = 1 if direction == 'forward' else 0
        args = args + ['--filterRNAstrand',direction,'--outFileFormat', 'bedgraph', '-bs', '1', '-p', '6','-o',bedgraphfilename,
                       '-r',whichchrom+":"+str(start)+'-'+str(end)]
        bamCoverage.main(args)
        file = open(bedgraphfilename,'rb')
        while True:
            line = file.readline()
            linelist = line.split()
            if (line == ''):
                break
            tmpstart = int(linelist[1])
            tmpend = int(linelist[2])
            data[tmpstart-start:tmpend-start,2*filenameidx+diridx] = float(linelist[3])
        file.close()

pyplot.bar([i-0.4 for i in range(data.shape[1])],sum(data))
pyplot.title(whichchrom+":"+str(start)+'-'+str(end))
pyplot.xticks(range(8),['AMF','AMR','APF','APR','RMF','RMR','RPF','RPR'])
