# Makes a TSS plot blindingly fast using MapReduce.
# For some reason this has to be run as a script, not "evaluated" to work.
# Similarly, *calls* to these functions cannot be made on the command line.
# Instead, you have to write some code that calls these functions and *run* them (from within PyCharm)
# Saving the output matrix to a file, "mapReduceHackeryOut.npy" which can then
# be plotting using the line at the bottom of this script.
#
# One problem with this is that reads that are falling off the edge of the window
# are not being counted.

import deeptools.mapReduce
import pysam
import numpy as np

# Tried to create a list of "tss" objects to pass around, but it
# didn't work because I was unable to pickle the list.
#class tss(object):
#    def __init__(self, a, b, c, d):
#        self.chrom = a
#        self.strand = b
#        self.pos = c
#        self.name = d

# For debugging
# fetch gives reads that extend beyond the boundaries given, which is good.
#counter = 0
#for read in bam.fetch('chr5',121641938,121641950) :
#   counter += 1
#   if read.flag in GOODFLAGS :
#       print(read.query_name, read.template_length, read.is_reverse, read.reference_start, read.next_reference_start, read.reference_end)

#       if (not read.is_reverse and read.template_length < 0) or (read.is_reverse and read.template_length > 0):
        #if read.query_name == 'HF5LMBBXX:6:2116:28452:12779':
#            print(read.query_name,read.template_length, read.is_reverse, read.reference_start)
#            print(read)
#            print('\n')
         #  break
#bam.close()

# IN PROGRESS!
# Input arguments to get_tss_alignedhist_fragment (need to be handed in as a tuple of 6 elements)
# 0) chrom
# 1) start
# 2) end
# 3) bam_file_name
# 4) maxreadlegnth
def get_read_length_hist_fragment(args):
    GOODFLAGS = {83, 99, 147, 163}
    readnameset = set()
    chrom, start, end, bam_file_name, maxreadlength = args
    bam = pysam.AlignmentFile(bam_file_name)

    profile = np.zeros(shape=(maxreadlength,1))
    for read in bam.fetch(chrom, start, end) : 
        if read.flag in GOODFLAGS : # Every read is represented twice in the BAM file
            if read.query_name in readnameset :
                readnameset.remove(read.query_name)
                continue
            else:
                readnameset.add(read.query_name)
                assert(read.is_proper_pair)
                len = abs(read.template_length)
                if len >= maxreadlength:
                    len = maxreadlength-1
                profile[len] += 1
    bam.close()
    return(profile)

# Histogram of read lengths    
def CalculateReadLengthHist (bamfilename, maxreadlength = 10000, chrom_sizes = None, blackListFileName = None):
    if chrom_sizes is None:
        bam = pysam.AlignmentFile(bamfilename, "rb")
        chrom_sizes = zip(bam.references, bam.lengths)
        bam.close()
    print ("calling mapreduce")
    result = deeptools.mapReduce.mapReduce((bamfilename, maxreadlength),
                                 get_read_length_hist_fragment,
                                 chrom_sizes,
                                 blackListFileName=blackListFileName, 
                                 genomeChunkLength=1000000,
                                 numberOfProcessors=6,
                                 verbose=True) # "result" is a list of ndarrays
    subtsshistograms =  np.concatenate(np.transpose(result))
    return(np.sum(subtsshistograms, axis = 1))


# Input arguments to get_tss_alignedhist_fragment (need to be handed in as a tuple of 6 elements)
# 0) chrom
# 1) start
# 2) end
# 3) bam_file_name
# 4) tsslist
# 5) halfwinwidth
def get_tss_alignedhist_fragment(args):
    GOODFLAGS = {83, 99, 147, 163}
    readnameset = set()
    chrom, start, end, bam_file_name, tsslist, halfwinwidth = args
    minitsslist = [t for t in tsslist if t[0] == chrom and t[2]>=start and t[2]<=end] # Only looking at TSSs between start and end
    bam = pysam.AlignmentFile(bam_file_name)

    profile = np.zeros(shape=(2*halfwinwidth,1))
    for tss in minitsslist:
        for read in bam.fetch(tss[0],tss[2]-halfwinwidth-100,tss[2]+halfwinwidth+100) :  # +/- 100 for edge effects?!
            if read.flag in GOODFLAGS : # Every read is represented twice in the BAM file
                if read.query_name in readnameset :
                    readnameset.remove(read.query_name)
                    continue
                else:
                    readnameset.add(read.query_name)
                if (tss[1]=='+'): # TSS is on the forward strand
                    if (not read.is_reverse):
                        a = read.reference_start - tss[2] + halfwinwidth
                        b = read.reference_start+abs(read.template_length) - tss[2] + halfwinwidth
                    else:
                        a = read.next_reference_start - tss[2] + halfwinwidth
                        b = read.reference_end - tss[2] + halfwinwidth
                else:
                    if (not read.is_reverse):
                        a = tss[2] - read.reference_start - abs(read.template_length) + halfwinwidth
                        b = tss[2] - read.reference_start + halfwinwidth
                    else:
                        b = tss[2] - read.next_reference_start + halfwinwidth
                        a = tss[2] - read.reference_end + halfwinwidth
                if (a>b):
                    print(a,b,read.is_reverse,read.template_length)
                    break
                assert(read.is_proper_pair)
                if (b < 0 or a > 2*halfwinwidth):
                    continue
                a = max(0,a)
                b = min(2*halfwinwidth,b)
                profile[a:b] += 1
    bam.close()
    return(profile)
 
# Histogram of reads aligned to TSSs    
def CalculateTSSHist(bamfilename, tsslist, halfwinwidth = 3000, chrom_sizes = None, blackListFileName = None):
    if chrom_sizes is None:
        bam = pysam.AlignmentFile(bamfilename, "rb")
        chrom_sizes = zip(bam.references, bam.lengths)
        bam.close()
    print ("calling mapreduce")
    result = deeptools.mapReduce.mapReduce((bamfilename, tsslist, halfwinwidth),
                                 get_tss_alignedhist_fragment,
                                 chrom_sizes,
                                 blackListFileName=blackListFileName, 
                                 genomeChunkLength=100000,
                                 numberOfProcessors=6,
                                 verbose=True) # "result" is a list of ndarrays
    subtsshistograms =  np.concatenate(np.transpose(result))
    return(np.sum(subtsshistograms, axis = 1))

def GetGapPos(filename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/gap.txt'):
    '''Function to make a GenomicInterval object that knows where all the gaps in the rhemac2.0 sequence.
    Returns a GenomicInterval.'''
    import HTSeq
    print "using gap file: {0}".format(filename)
    gapfile = open(filename)
    gappos = HTSeq.GenomicArray("auto", stranded=False)
    while True:
        line = gapfile.readline()
        linelist = line.split()
        if line == '':
            break
        c = HTSeq.GenomicInterval(linelist[0], int(linelist[1]), int(linelist[2]))
        gappos[c] = 1
    gapfile.close()
    return(gappos)

def EliminateGapsFromTssList(tsslist, halfwinwidth, filename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/gap.txt'):
    # filename (of a bed file of gaps) is optional
    import HTSeq
    newtsslist = list()
    gappos = GetGapPos(filename)
    for i in tsslist:
        if (i[2] - halfwinwidth > 0): # Should also have a check for going beyond the end of the chromosome
            c = HTSeq.GenomicInterval(i[0], i[2]-halfwinwidth, i[2]+halfwinwidth)
            steps = [j for j in gappos[c].steps()]
            if len(steps) == 1:
                newtsslist.append(i)
    return newtsslist

# Counting reads from the BAM file that overlap with NUMTs
def CountBAMOverlapReadsHelper(args):
    GOODFLAGS = {83, 99, 147, 163}
    chrom, start, end, bam_file_name, regionlist = args
    miniregionlist = []
    for t in regionlist:
        if t[0] == chrom:
            if (t[2]>=start and t[3]<=end): # regions that begin and end between 'start' and 'stop'
                miniregionlist.append(t)
            elif (t[2] < start and t[3] > end): # regions that begin before 'start' and end after 'stop'
                miniregionlist.append([t[0],t[1],start,end])
            elif (t[2] < start and t[3] <= end): # regions that begin before 'start' and end before 'stop'
                miniregionlist.append([t[0],t[1],start,t[3]])
            elif (t[2] < start and t[3] <= end):  # regions that begin after 'start' and end after 'stop'
                miniregionlist.append([t[0], t[1], t[2], end])
    bam = pysam.AlignmentFile(bam_file_name)
    counter = 0
    for region in miniregionlist:
        for read in bam.fetch(region[0],region[2],region[3]) :
            if read.flag in GOODFLAGS : # Every read is represented twice in the BAM file
                counter += 1
    bam.close();            
    return counter

# Counting the number of BAM reads that overlap with region list.
def CountBAMOverlapReads(bamfilename, regionlist, chrom_sizes = None):
    if chrom_sizes is None:
        bam = pysam.AlignmentFile(bamfilename, "rb")
        chrom_sizes = zip(bam.references, bam.lengths)
        bam.close()
    result = deeptools.mapReduce.mapReduce((bamfilename, regionlist),
                                 CountBAMOverlapReadsHelper,
                                 chrom_sizes,
                                 genomeChunkLength=10000000,
                                 numberOfProcessors=6,
                                 verbose=True) # "result" is a list of ndarrays
    return(np.sum(result))

# Find a data file in the directory tree
def findfile(name, path = '/Volumes/2TBdrive'):
    import os
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)