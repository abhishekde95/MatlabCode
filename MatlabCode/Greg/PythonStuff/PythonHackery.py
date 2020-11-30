# This is my flailing around with the data structures and functions in HTSeq.py.
# There's a function, tss_aligned_hist, that makes a histogram aligned to annotated
# transcription start sites.
# There's a function, MakeRandomTSSGenomicArrayOfSets, that randomizes the positions of
# the transcription start sites within each chromosome (randomizing the inter-TSS-intervals)
# and returns them as a GenomicArrayOfSets, which is one of HTSeqs useful data structures.
# There's a bit of code at the end of this file that obtains the number of ATAC-seq reads
# in the neighborhood of TSSs in two different conditions (e.g. M and P) and compares
# these to the number of counts from RNA-seq datat file.
#
# Finally making this into a function so I can iterate over
# things and call this function again and again.
# "lengths" only looks at fragments that were within 2*halfwinwidth of a TSS.

class BAMFileSortError(Exception):
    pass

# Makes a histogram of reads aligned to the TSS where + numbers mean downstream (3') and - numbers mean upstream (5').
# bamfile has to be sorted by read name, otherwise HTSeq.pair_SAM_alignments crashes.
# tsspos can either be a GenomicArrayofSets or a list of GenomicArrayofSets. The latter will be
# useful for making multiple histograms with only a single pass through the BAM file.
def tss_aligned_hist (bamfile, tsspos, halfwinwidth = 3000, sortreadsbysize=True, makefigure = True, howmanyreads = None):
    # hitting a buffer overloaded error unless I use files that are sorted by name

    GOODFLAGS = {83, 99, 147, 163}
    LONGREAD = 250  # bp that is dividing line between long and medium reads
    SHORTREAD = 100  # bp that is dividing line between medium and short reads

    # making tsspos a list, even if it only consists of one element
    if type(tsspos) == HTSeq.GenomicArrayOfSets:
        tsspos = [tsspos,]

    read_seq_iter = iter(bamfile)
    durablebamreader = HTSeq.pair_SAM_alignments(read_seq_iter)  # Should have been doing it this way all along
    profile = np.zeros( (2*halfwinwidth,3,len(tsspos)), dtype="i" )
    counter = 0 # counting PE reads from BAM file
    lengths = dict()
    print('In TssAlignedHist')
    try:
        if howmanyreads is None:
            myiterator = durablebamreader
        else:
            myiterator = itertools.islice(durablebamreader,howmanyreads) # For debugging
        for a, b in myiterator:
            counter += 1
            if a is None or b is None:
                #raise BAMFileSortError
                if (a is None):
                    a = b
                    a.flag = -1
                if (b is None):
                    b = a
                    b.flag = -1
            elif a.flag in GOODFLAGS and b.flag in GOODFLAGS:
                s = set()
                if (a.iv.strand == '+'):
                    c = HTSeq.GenomicInterval(a.iv.chrom, a.iv.start, b.iv.end) # spanning reads a and b
                    alignment_length = b.iv.end - a.iv.start
                    if a.iv.start > b.iv.start :
                        print('start and ends in wrong order? (1)')
                    if alignment_length < 0 :
                        print ('alignment length < 0! (1)')
                       # break
                elif (a.iv.strand == '-'): # b.iv.strand = '+'
                    c = HTSeq.GenomicInterval(a.iv.chrom, b.iv.start, a.iv.end)  # spanning reads a and b
                    alignment_length = a.iv.end-b.iv.start
                    if b.iv.start > a.iv.start :
                        print('start and ends in wrong order? (2)')
                    if alignment_length < 0:
                        print ('alignment length < 0! (2)')
                              # break
                for whichtsspos in range(len(tsspos)):
                    tmptsspos = tsspos[whichtsspos]
                    s = set()
                    for step_iv, step_set in tmptsspos[ c ].steps():
                        s |= step_set # step_set is set of GenomicPosition objects, each with with length 1
                    for p in s:
                        if p.strand == "+": # c is set up as if it's on the '+' strand
                            start_in_window = c.start - p.pos + halfwinwidth  # This "+halfwinwidth" takes -3000 (which is legit) and moves it to index 0
                            end_in_window   = c.end   - p.pos + halfwinwidth
                        else:  # Reversing the orientation of alignments on the - strand so that
                            #  [500 510] when p.pos = 508 becomes [-2 8]. Then we add halfwinwidth so that 0 is the first index.
                            start_in_window = p.pos + halfwinwidth - c.end
                            end_in_window   = p.pos + halfwinwidth - c.start
                        start_in_window = max( start_in_window, 0 )
                        end_in_window = min( end_in_window, 2*halfwinwidth )
                        if (alignment_length <= SHORTREAD):
                            profile[ start_in_window : end_in_window, 0, whichtsspos] += 1
                        elif (alignment_length <= LONGREAD):
                            profile[ start_in_window : end_in_window, 1, whichtsspos] += 1
                        else:
                            profile[ start_in_window : end_in_window, 2, whichtsspos] += 1
                if alignment_length in lengths: # keeping track of alignment lengths in a dictionary
                    lengths[alignment_length] += 1
                else:
                    lengths[alignment_length] = 1
            if counter % 10000 == 0:
                print (counter)
    except ValueError as inst :
        print("Hit an error")
        print(type(inst))
        print(inst.args)
        print(inst)
        print(a)
        print(b)
        print(bamfile.get_line_number_string())
    except BAMFileSortError :
        print ('At least one of two PE reads is "None"')
        print ("You're probably using a coordinate-sorted BAM file instead of a name-sorted one")
        import pdb
        pdb.set_trace()
    except KeyboardInterrupt :
        print('Ctrl-C detected. Type "exit" to return to exit the debugger.')
        import pdb
        pdb.set_trace()
    except :
        print('Unknown error')
        errinfo = sys.exc_info()
        print(errinfo)
        import pdb
        pdb.set_trace()

    if makefigure: # Pretty outdated
        pyplot.figure()
        pyplot.subplot(2,1,1)
        pyplot.title(bamfile.filename[bamfile.filename.rfind(os.path.sep) + 1:])
        if (sortreadsbysize):
            pyplot.plot(np.arange( -halfwinwidth, halfwinwidth ), profile)
            pyplot.legend(['<' + str(SHORTREAD), '{0}-{1}'.format(str(LONGREAD), str(SHORTREAD)), '>' + str(LONGREAD)], loc = 2, fontsize='small')
        else:
            pyplot.plot(np.arange(-halfwinwidth, halfwinwidth), np.sum(profile, axis = 1))
        ax = pyplot.subplot(2,1,2)
        pyplot.plot(lengths.keys(), lengths.values(),'-')
        if sortreadsbysize:
            pyplot.plot([SHORTREAD, SHORTREAD],[10, 10000],'b-')
            pyplot.plot([LONGREAD, LONGREAD],[10, 10000],'g-')
        ax.set_xscale('log')
        ax.set_yscale('log')
        pyplot.ylim(10, 10000)
        ax.set_ylabel('count')
        ax.set_xlabel('fragment length (bp)') # plot color order: blue, green, red
    return(profile)

##############################################
# Randomizing the transcription start sites on each chromosome
def MakeRandomTSSGenomicArrayOfSets(tsspos_set, halfwinwidth = 3000):
    from random import shuffle
    uniquechroms_chr = set([])
    for p in tsspos_set:
        if not p.chrom[0].isdigit():
            uniquechroms_chr.add(p.chrom)

    chromdict = {} # dict
    for p in tsspos_set:
        if p.chrom in uniquechroms_chr: # including chrMT for now
            if p.chrom in chromdict:
                chromdict[p.chrom].append(p.start)
            else:
                chromdict[p.chrom] = [p.start]
    for k in chromdict.keys():
        chromdict[k].sort() # in place sort
        ItssIs = list(np.diff([1,]+chromdict[k]))
        shuffle(ItssIs)  # in place
        chromdict[k] = list(np.cumsum(ItssIs))

    # now have to make a GenomicArrayofSets with (potentially overlapping) windows of length 2*halfwinwidth
    tsspos_fake = HTSeq.GenomicArrayOfSets( "auto", stranded=False ) # randomized TSSs
    for p in sorted(tsspos_set, key=lambda p: p.chrom) :
        try:
            p.pos = chromdict[p.chrom].pop(0)
            if (p.pos > halfwinwidth):
                window = HTSeq.GenomicInterval(p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, ".")
                tsspos_fake[window] += p
        except KeyError:
            pass
        except:
            print("unexpected error")
            raise
    return (tsspos_fake)

def GetGapPos(filename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/gap.txt'):
    '''Function to make a GenomicInterval object that knows where all the gaps in the rhemac2.0 sequence.
    Returns a GenomicInterval.'''
    import HTSeq
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
###########################################################################
# sorting genes into expressed and non-expressed. Then looking for ATAC-seq
# peaks at the TSSs. Don't forget, we need to express gene count in RPKM
# or something like that.
# 9/25/16 updating the plot on every iteration isn't working well.
# Making tss_aligned_hist() return a matrix each column of which is a
# tss-aligned histogram. Then I can plot them all at once from the calling
# function which will hopefully crash less.
###########################################################################
# My hacking around with Python
import sys
import HTSeq
from matplotlib import pyplot
import numpy as np
import itertools
import os
import pickle

# Here's the command you need to get this file in the python path
# Except that it didn't seem to work
#sys.path + ['/Users/horwitzlab/Desktop/SEQ_ANALYSIS/customtools']
bamfilename = '/Volumes/2TBdrive/ATAC/ATAC_2016_12/Z13145_V1t_q_namesorted.bam'
#bamfilename = '/Volumes/2TBdrive/ATAC/ATAC_2016_09/Z12131_M_rhemac8.bam'
RANDOMIZE_TSS_POSITIONS = False # outdated
SKIPGENESNEARGAPS = False
bamfile = HTSeq.BAM_Reader(bamfilename) # It would be great to figure out how to set the path that BAM_Reader will check for files.
gtffile = HTSeq.GFF_Reader("/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf")

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

genehitprctiles = [50, 60, 70, 80, 90]
genehitprctiles = [0,]
cutoffvalues = np.percentile(RNAdict.values(),genehitprctiles)
cutoffvalues = np.append(cutoffvalues, maxreads)

# Getting the ranks of the RNA expression levels (using raw read counts)
RNAranks = {}
for a, b in RNAdict.iteritems():
    for j in range(len(cutoffvalues)):
       if a not in RNAranks and b <= cutoffvalues[j] :
           RNAranks[a] = j # ranks of gene expression
           break
#pyplot.figure()
#pyplot.hist(RNAranks.values())
#pyplot.show()

halfwinwidth = 3000
if SKIPGENESNEARGAPS :
    gappos = GetGapPos()
else :
    gappos = None
tssposlist = [] # A list of GenomicArrayOfSets objects - one for each level of RNA expression used to separate genes
print('Getting TSS positions for each level of expression')
for expressionrankcounter in range(len(cutoffvalues)): # This would go faster if we weren't making n passes through the GTF file
    print('expressionrankcounter is {0}'.format(expressionrankcounter))
    tsspos = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    # if (RANDOMIZE_TSS_POSITIONS): # Doesn't deal with gaps in the sequence!
    #     tsspos_set = set() # tsspos_set is a set, tsspos is a GenomicArrayOfSets
    #     for feature in gtffile:
    #         if feature.type == "exon" and feature.attr["exon_number"] == "1" and RNAranks[feature.name] == expressionrankcounter :
    #             tsspos_set.add(feature.iv.start_d_as_pos)
    #     tsspos = MakeRandomTSSGenomicArrayOfSets(tsspos_set, halfwinwidth=halfwinwidth)
    # else:
    for feature in gtffile:
        if feature.type == "exon" and feature.attr["exon_number"] == "1" and RNAranks[feature.name] == expressionrankcounter :
            p = feature.iv.start_d_as_pos # start_d_as_pos = The position of the TSS irrespective of strand
            if (p.start > halfwinwidth):
                window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
                if SKIPGENESNEARGAPS :
                    steps = [i for i in gappos[window].steps()]
                    if len(steps) > 1:  # There should only be one "step" in c if there's no gap near the TSS
                        continue
                tsspos[ window ] += p # Cool. Each window of tsspos is of length 2*halfwinwidth
    tssposlist.append(tsspos)
print('Calling TssAlignedHist')
tsshistmat = tss_aligned_hist (bamfile, tssposlist, halfwinwidth, sortreadsbysize=True, makefigure = False, howmanyreads = None) # Include howmanyreads = 50000 for debugging
pickle.dump(tsshistmat,open('tsshistmat.p','wb'))

# Summing across read lengths
pyplot.figure()
ax = pyplot.axes()
for i in range(tsshistmat.shape[2]):
    pyplot.plot(np.arange(-halfwinwidth, halfwinwidth),np.sum(tsshistmat[:,:,i], axis = 1))
filenamestr = bamfile.filename[bamfile.filename.rfind(os.path.sep) + 1:]
pyplot.title(filenamestr)
#
legendstr = ['< '+str(genehitprctiles[0])]
if len(genehitprctiles) > 2:
   legendstr += [str(genehitprctiles[i])+'-'+str(genehitprctiles[i+1])+'%' for i in range(len(genehitprctiles)-1)]
legendstr += ['> '+str(genehitprctiles[-1])]
pyplot.legend(legendstr,fontsize = 10)

# Summing across expression levels
pyplot.figure()
ax = pyplot.axes()
for i in range(tsshistmat.shape[1]):
    pyplot.plot(np.arange(-halfwinwidth, halfwinwidth),np.sum(tsshistmat[:,i,:], axis =1 ))
filenamestr = bamfile.filename[bamfile.filename.rfind(os.path.sep) + 1:]
pyplot.title(filenamestr+ ' by read length')

# if RANDOMIZE_TSS_POSITIONS:
#     filenamestr += '.RandomTSSs'
pyplot.savefig(filenamestr + str(expressionrankcounter)+'.pdf')


####################################
# Getting counts of ATAC-seq reads within 1000 bp or so of TSS for both
# magno and parvo. Then comparing this to RNA-seq hits.
####################################
import sys
import HTSeq
from matplotlib import pyplot
import numpy as np
import itertools
import os
import pickle
ATACseqpath = '/Volumes/2rive/ATAC/ATAC_2016_09/'
ATACseqfilenames = ['Z12131_M_namesorted.bam','Z12131_P_namesorted.bam']
RNAseqpath = '/Volumes/2TBdrive/RNA/RNA_2016_01'
RNAseqfilenames = ['A09177_M_q_genecounts','A09177_P_q_genecounts']
RNAseqpath = '/Volumes/2TBdrive/RNA/RNA_2016_05'
RNAseqfilenames = ['A20416_M_q_genecounts','A20416_P2_q_genecounts'] # Add counts from other RNA-seq runs?

GOODFLAGS = {83, 99, 147, 163}
SAVEWORKSPACE = True # Saving tsspos, RNAdict, and ATACdict in case you hit a crash
halfwinwidth = 1000

# Making a dictionary of gene names as keys and lists as values where the first element in each
# list is the count from the first RNAseqfile and the second element is the count from the second
RNAdictlist = []
for filename in RNAseqfilenames:
    if not os.path.isfile(RNAseqpath+os.path.sep+filename):
        my_error = IOError("{0} does not seem to exist.".format(RNAseqpath+os.path.sep+filename))
        raise my_error
    file = open(RNAseqpath+os.path.sep+filename)
    RNAdict={}
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
    RNAdictlist.append(RNAdict)
# Code below assumes that there are only two files
for k in RNAdictlist[0].keys(): RNAdictlist[0][k] = (RNAdictlist[0][k], RNAdictlist[1][k])
RNAdict = RNAdictlist[0]
del RNAdictlist

# Can use max counts to normalize gene hits
maxcounts = np.zeros(shape = (2,1),dtype = 'int')
for k in RNAdict:
    maxcounts = maxcounts+np.asarray(RNAdict[k])[:,None]

# Trying plotting
# from __future__ import division
# import time
#
# pyplot.ioff()
# pyplot.figure()
# start_time = time.time()
# # This calculation, below takes a full minute
# for k in RNAdict:
#     pyplot.plot(RNAdict[k][0]/maxcounts[0],RNAdict[k][1]/maxcounts[1],'k.')
# end_time = time.time()
# print('elapsed time is {0}'.format(end_time-start_time))
# start_time = time.time()
# pyplot.show() # This takes forever - don't do it!
# end_time = time.time()
# ax = pyplot.gca()
# ax.set_yscale('log')
# ax.set_xscale('log')

# Now setting up the tsspos GenomicArray so that we can run through and see how many ATAC-seq
# fragments fall within ±1 kb or so.
gtffile = HTSeq.GFF_Reader("/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf")
tsspos = HTSeq.GenomicArrayOfSets("auto", stranded=False)
for feature in gtffile:
    if feature.type == "exon" and feature.attr["exon_number"] == "1" :
        p = feature.iv.start_d_as_pos  # start_d_as_pos = The position of the TSS irrespective of strand
        if (p.start > halfwinwidth):
            window = HTSeq.GenomicInterval(p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, feature.iv.strand)
            tsspos[window] += (feature.name, p.pos)  # replacing 'p' (a GenomicPosition) with a tuple consisting of
            # the feature name (a string) followed by the position in the genome. This way I can keep track of
            # both the position of the TSS and the name of the gene.

# Now iterating over the reads to find one that overlaps with a TSS
ATACdict = {key: np.zeros((2,1)) for key in RNAdict.keys()} # Initializing a dict with the same keys as the RNAdict
print('Counting ATAC-seq reads')
alldone = False # for debugging
for filename in ATACseqfilenames:
    counter = 0
    bamfile = HTSeq.BAM_Reader(ATACseqpath+filename)
    read_seq_iter = iter(bamfile)  # this does work
    first_read = read_seq_iter.next()
    read_seq = itertools.chain([first_read], read_seq_iter)
    durablebamreader = HTSeq.pair_SAM_alignments(read_seq)
    for a, b in durablebamreader:
    #for a, b in itertools.islice(durablebamreader,900000):  # for debugging
        counter += 1
        if a is None or b is None:
            print("Warning! Discarding a read that does not appear to have a mate.")
            import pdb
            pdb.set_trace()
        elif a.flag in GOODFLAGS and b.flag in GOODFLAGS:
            s = set()
            if (a.iv.strand == '+'):
                c = HTSeq.GenomicInterval(a.iv.chrom, a.iv.start, b.iv.end)  # spanning reads a and b
            elif (a.iv.strand == '-'):  # b.iv.strand = '+'
                c = HTSeq.GenomicInterval(a.iv.chrom, b.iv.start, a.iv.end)  # spanning reads a and b
            for step_iv, step_set in tsspos[c].steps():
                s |= step_set
            for p in s:
                ATACdict[p[0]][ATACseqfilenames.index(filename)] += 1
            alldone = False # For debugging
        if counter % 10000 == 0:
            print (counter)
    if alldone: # For debugging
        break

if SAVEWORKSPACE:
    workspace = [tsspos, RNAdict, ATACdict]
    pickle.dump(workspace,open('workspace.p','wb'))
if LOADWORKSPACE:
    workspace = pickle.load(open("workspace.p", 'rb'))  # In case you hit a crash
    tsspos = workspace[0]
    RNAdict = workspace[1]
    ATACdict = workspace[2]

# Not at all clear that RNAarray and ATACarray are in the same order
RNAarray = np.squeeze(np.asarray([RNAdict[k] for k in ATACdict.keys()]))
RNAarray = np.squeeze(np.asarray(RNAdict.values()))

ATACarray = np.squeeze(np.asarray([ATACdict[k] for k in ATACdict.keys()]))
pyplot.figure()
pyplot.plot(ATACarray[:,0],ATACarray[:,1],'.')
ax = pyplot.gca()
ax.set_yscale('log')
ax.set_xscale('log')

from __future__ import division
[i for i, j in ATACdict.items() if sum(j) > 0 and abs(j[0]-j[1])/(j[0]+j[1]) > .5 and np.max(j)>100] # Differential ATAC-seq peaks


# This crashes!
pyplot.plot((ATACarray[:,0]-ATACarray[:,1])/(ATACarray[:,0]+ATACarray[:,1]),
            (RNAarray[:,0] - RNAarray[:,1])/(RNAarray[:, 0] + RNAarray[:, 1]),'.')
pyplot.xlabel('ATAC-seq (M-P)/(M+P)')
pyplot.ylabel('RNA-seq (M-P)/(M+P)')

pyplot.plot(ATACarray[:,0]-ATACarray[:,1], RNAarray[:,0] - RNAarray[:,1],'.')

np.corrcoef(np.append(ATACarray, RNAarray,axis = 1),rowvar = 0)
np.corrcoef(ATACarray[:,0]-ATACarray[:,1], RNAarray[:,0]-RNAarray[:,1])

[i for i, j in ATACdict.items() if abs(j[0]-j[1]) > 100 and abs(RNAdict[i][0]-RNAdict[i][1]) > 1000]

# More ATAC-seq peaks piling up near TSSs for magno sample?
from scipy.stats import ttest_1samp
t, p = ttest_1samp(ATACarray[:,0]-ATACarray[:,1],0)
np.mean(ATACarray,axis=0)
#########################################################################################
#
# Loading in a bunch of RNA-seq files and calculate something like a F-statsitic
# on them.
#
# #########################################################################################

RNAseqpathroot = '/Volumes/2TBdrive/RNA'
RNAseqsubpaths = ['RNA_2016_01', 'RNA_2016_05']
RNAseqfilenames = ['A09177_M_q_genecounts', 'A09177_P_q_genecounts',
                   'A20416_P2_q_genecounts', 'A21015_M_q_genecounts',
                   'A20416_M_q_genecounts','A20416_P1_q_genecounts']
Magnofiles = {'A09177_M_q_genecounts','A21015_M_q_genecounts','A20416_M_q_genecounts'}
Parvofiles = {'A09177_P_q_genecounts','A20416_P2_q_genecounts','A20416_P1_q_genecounts'}
Malefiles  = {'A21015_M_q_genecounts','A09177_M_q_genecounts','A09177_P_q_genecounts'}
Femalefiles = {'A20416_M_q_genecounts','A20416_P1_q_genecounts','A20416_P2_q_genecounts'}

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
    RNAdictlist.append(tmpRNAdict)

RNAdict = {}
for key in RNAdictlist[0].keys():
    tmplist = []
    for i in range(len(RNAdictlist)):
        tmplist.append(RNAdictlist[i][key])
        del RNAdictlist[i][key]
    RNAdict[key] = tmplist[:]  # deep copy since I keep manipulating tmplist
for i in RNAdictlist: assert not i, "Extra genes detected in one or more RNA dictionaries."
RNAarray = np.array(RNAdict.values(), dtype = int)
geneNames = RNAdict.keys()
# Getting column indices M vs P and Male vs Female
magnoidxs = [i for i, j in enumerate(RNAseqfilenames) if j in Magnofiles]
parvoidxs = [i for i, j in enumerate(RNAseqfilenames) if j in Parvofiles]
maleidxs = [i for i, j in enumerate(RNAseqfilenames) if j in Malefiles]
femaleidxs = [i for i, j in enumerate(RNAseqfilenames) if j in Femalefiles]
assert not set(magnoidxs).intersection(parvoidxs), "Some magno are also parvo"
assert not set(maleidxs).intersection(femaleidxs), "Some males are also female"

# normalizing to proportion of total reads per sample
normRNAarray1 = RNAarray.astype(float)/np.sum(RNAarray,axis = 0)[None,:]
normRNAarray2 = np.sqrt(RNAarray)/np.sum(np.sqrt(RNAarray),axis = 0)[None,:]
normRNAarray = normRNAarray2 # Allows user to adjust the normalization we're using
# Should we divide by gene lengths?
# One problem with this is that we're weighting the MiSeq data the same as the HiSeq data
F = np.zeros(normRNAarray1.shape[0])
ns = np.zeros(normRNAarray1.shape[0])
dumbstat = np.zeros(normRNAarray1.shape[0]) # Just the difference in the total number of (normalized) reads
COMPAREGENDERS = True
if COMPAREGENDERS:
    grp0idxs, grp1idxs = femaleidxs, maleidxs
else:
    grp0idxs,grp1idxs = magnoidxs, parvoidxs

for i in range(normRNAarray.shape[0]):
    grp0 = normRNAarray[i][grp0idxs]
    grp1 = normRNAarray[i][grp1idxs]
    withinGroupVar = (sum((grp0-np.mean(grp0))**2)+
                      sum((grp1-np.mean(grp1))**2))/(normRNAarray.shape[1]-2)
    betweenGroupVar = (grp0.shape[0]*(np.mean(grp0)-np.mean(normRNAarray[i]))**2 +
                       grp1.shape[0]*(np.mean(grp1)-np.mean(normRNAarray[i]))**2)
    F[i] = betweenGroupVar/withinGroupVar
    ns[i] = sum(normRNAarray[i])
    #dumbstat[i] = (np.mean(grp0)-np.mean(grp1))/(np.mean(grp0)+np.mean(grp1))
#F = dumbstat # dumbstat seems to work pretty well
fig1 = pyplot.figure()
pyplot.plot(F,ns,'.',markersize=2)
pyplot.gca().set_xscale('log')
pyplot.gca().set_yscale('log')
pyplot.xlabel('F')
pyplot.ylabel('n')

geneidxs = np.logical_and(F>np.percentile(F[~np.isnan(F)],99), ns > np.median(ns))
from itertools import compress
interestingCandidates = list(compress(geneNames, geneidxs))
print("number of interesting candidates: {0}".format(sum(geneidxs)))

if not len(interestingCandidates) == 0:
    for key in interestingCandidates:
        genecounts = RNAdict[key]
        #print (key, [normRNAarray1[i] for i in grp0idxs], [normRNAarray1[i] for i in grp1idxs])
        print (key, [genecounts[i] for i in grp0idxs], [genecounts[i] for i in grp1idxs])

# Looking at correlations across the RNA samples
fig2 = pyplot.figure()
S2 = np.corrcoef(RNAarray[:,femaleidxs+maleidxs], rowvar = 0)
ax = fig2.add_subplot(2,2,1)
pyplot.imshow(S2, interpolation = 'none')
ax.set_title('Female vs Male')

S2 = np.corrcoef(RNAarray[:,magnoidxs+parvoidxs], rowvar = 0)
ax = fig2.add_subplot(2,2,2)
pyplot.imshow(S2, interpolation = 'none')
ax.set_title('Magno vs Parvo')

# Looking at the correlations among the three statistics I've been playing with
# (dumbstatistic, F (normalized to total read count), F (sqrt counts and then normalized to total)


#########################################################################################
#
# Reading in a text file made macs2 and looks for ATAC-seq peaks that
# occur withing 2 kb or so or a TSS. At some point I should also look for
# ATAC-seq peaks within a few bp of gaps, because that seems to happen very
# often.
#
#########################################################################################

# First, setting up the genomic intervals
halfwinwidth = 2000
gtffile = HTSeq.GFF_Reader("/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf")
tsspos = HTSeq.GenomicArrayOfSets("auto", stranded=False)
for feature in gtffile:
    if feature.type == "exon" and feature.attr["exon_number"] == "1" :
        p = feature.iv.start_d_as_pos  # start_d_as_pos = The position of the TSS irrespective of strand
        if __name__ == '__main__':
            if (p.start > halfwinwidth):
                window = HTSeq.GenomicInterval(p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, feature.iv.strand)
                tsspos[window] += (feature.name, p.pos)  # replacing 'p' (a GenomicPosition) with a tuple consisting of
                # the feature name (a string) followed by the position in the genome. This way I can keep track of
                # both the position of the TSS and the name of the gene.

# I think  the "narrowPeak" file is the one I want. Here's the contents by column
# 1) chromosome (e.g. "chr1")
# 2) start position
# 3) end position
# 4) name of the peak (e.g. "/Users/horwitzlab/Desktop/SEQ_DATA/ATAC/ATAC_2016_05/Z12125_P2_q_peaks_peak_1")
# 5) "integer score for display"? A UCSC genome browser specific thing
# 6) "."
# 7) fold-change (relative to some kind of baseline)
# 8) -log10(pvalue)
# 9) -log10(qvalue)
# 10) Position of summit relative to peak start
MACS2outputfilename = 'Z12131_M_q_peaks_peaks.narrowPeak'
MACS2outputpath = '/Volumes/2TBdrive/ATAC/ATAC_2016_09'
interesting_genes = []
try:
    file = open(MACS2outputpath+os.path.sep+MACS2outputfilename)
except IOError as inst:
    print(MACS2outputpathpath+os.path.sep+MACS2outputfilename+' does not seem to exist')
else:
    while True:
        line = file.readline()
        if line == '':
            break
        linelist = line.split()
        c = HTSeq.GenomicInterval(linelist[0], int(linelist[1]), int(linelist[2]))
        s = set()
        for step_iv, step_set in tsspos[c].steps():
            s |= step_set
        if len(s) > 0:
            tmplist = []
            print(linelist[0]+':'+linelist[1]+'-'+linelist[2])
            tmplist.append(list(s)[0][0]) # gene name
            tmplist.append(int(linelist[1])) # start
            tmplist.append(int(linelist[2]))  # stop
            tmplist.append(float(linelist[7]))  # fold-change
            tmplist.append(float(linelist[8]))  # -log10 p-value
            tmplist.append(float(linelist[9]))  # -log10 q-value
            interesting_genes.append(tmplist)
    file.close()
# Comparing read counts to ATAC-seq peaks
ATACpeakgenenames = set()
for l in interesting_genes:
    ATACpeakgenenames.add(l[0])

RNAhits_near_ATAC_peaks =list()
RNAhits_farfrom_ATAC_peaks =list()

for k in RNAdict.keys():
    if (k in ATACpeakgenenames):
        print('got one')
        RNAhits_near_ATAC_peaks.append(np.mean(RNAdict[k]))
    else:
        RNAhits_farfrom_ATAC_peaks.append(np.mean(RNAdict[k]))

print(np.mean(RNAhits_near_ATAC_peaks), np.mean(RNAhits_farfrom_ATAC_peaks))
print(np.std(RNAhits_near_ATAC_peaks), np.std(RNAhits_farfrom_ATAC_peaks))
print(len(RNAhits_near_ATAC_peaks), len(RNAhits_farfrom_ATAC_peaks))
from scipy.stats import ttest_ind
t, p = ttest_ind(RNAhits_near_ATAC_peaks, RNAhits_farfrom_ATAC_peaks)



#########################################################################################
#
# Counting the number of reads that overlap some fixed window that is specified relative
# to the TSS. Can I do this with both deeptools HTseq libraries and get the same answer?
# Need to skip TSSs near gaps!
#
# Hopefully I can figure out a way to leverage deeptools "mapReduce" function to make
# this go fast.
#
#########################################################################################
# Getting gaps
gapfilename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/gap.txt'
gappos = GetGapPos(gapfilename)

# First using HTseq
gtffile = HTSeq.GFF_Reader("/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf")
halfwinwidth = 1000
tsspos = HTSeq.GenomicArrayOfSets("auto", stranded=False)
onekb_upstreampos = HTSeq.GenomicArrayOfSets("auto", stranded=False)
twokb_upstreampos = HTSeq.GenomicArrayOfSets("auto", stranded=False)
ATACdict = dict()

for feature in gtffile:
    if feature.type == "exon" and feature.attr["exon_number"] == "1":
        p = feature.iv.start_d_as_pos  # start_d_as_pos = The position of the TSS irrespective of strand
        c = HTSeq.GenomicInterval(p.chrom, p.pos - 2000 - halfwinwidth, p.pos + halfwinwidth,'.')
        if c.start < 0:
            continue
        steps = [i for i in gappos[c].steps()]
        if len(steps) > 1: # There should only be one "step" in c is there's no gap near the TSS
            #print("skipping {} because number of steps is {}".format(feature.name, len(steps)))
            continue
        else:
            #print("keeping {}".format(feature.name))
        ATACdict[feature.name] = np.zeros(shape = (3,1))
        if (p.start > halfwinwidth): # Just to avoid TSSs that are at the very begining of chromosomes
            window = HTSeq.GenomicInterval(p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, ".")
            tsspos[window] += (feature.name, p.pos)
        if (p.start > halfwinwidth + 1000):
            window = HTSeq.GenomicInterval(p.chrom, p.pos - 1000 - halfwinwidth, p.pos - 1000 + halfwinwidth, ".")
            onekb_upstreampos[window] += (feature.name, p.pos)
        if (p.start > halfwinwidth + 2000):
            window = HTSeq.GenomicInterval(p.chrom, p.pos - 2000 - halfwinwidth, p.pos - 2000 + halfwinwidth, ".")
            twokb_upstreampos[window] += (feature.name, p.pos)

bamfilename = '/Volumes/2TBdrive/ATAC/ATAC_2016_09/Z12131_V1b_namesorted.bam'
bamfile = HTSeq.BAM_Reader(bamfilename)
read_seq_iter = iter(bamfile)
durablebamreader = HTSeq.pair_SAM_alignments(read_seq_iter)

GOODFLAGS = {83, 99, 147, 163}
howmanyreads = None
counter = 0
if howmanyreads is None:
    myiterator = durablebamreader
else:
    myiterator = itertools.islice(durablebamreader, howmanyreads)  # For debugging
try:
    for a, b in myiterator:
        counter += 1
        if a is None or b is None:
            raise BAMFileSortError
        elif a.flag in GOODFLAGS and b.flag in GOODFLAGS:
            if (a.iv.strand == '+'):
                c = HTSeq.GenomicInterval(a.iv.chrom, a.iv.start, b.iv.end)  # spanning reads a and b
            elif (a.iv.strand == '-'):  # b.iv.strand = '+'
                c = HTSeq.GenomicInterval(a.iv.chrom, b.iv.start, a.iv.end)  # spanning reads a and b

            s = set()
            for step_iv, step_set in tsspos[c].steps():
                s |= step_set
            for p in s:
                ATACdict[p[0]][0] += 1

            s = set()
            for step_iv, step_set in onekb_upstreampos[c].steps():
                s |= step_set
            for p in s:
                ATACdict[p[0]][1] += 1

            s = set()
            for step_iv, step_set in twokb_upstreampos[c].steps():
                s |= step_set
            for p in s:
                ATACdict[p[0]][2] += 1

        if counter % 10000 == 0:
            print (counter)
except KeyboardInterrupt:
    raise
workspace = [ATACdict]
pickle.dump(workspace, open('workspace.p', 'wb'))
workspace = pickle.load(open("workspace.p", 'rb'))  # In case you hit a crash
ATACdict = workspace[0]

readcountarray = np.zeros(shape = (len(ATACdict),3))
counter = 0
for k in ATACdict.keys():
    readcountarray[counter,:] = np.squeeze(ATACdict[k])
    counter += 1

np.corrcoef(readcountarray,rowvar=0)
pyplot.plot(readcountarray[:,0],readcountarray[:,2],'.')
pyplot.show()

for k in ATACdict.keys():
    if __name__ == '__main__':
        if ATACdict[k][0] >25 and ATACdict[k][1] < 10  and ATACdict[k][2] > 25:
            print(k)
         
#%%
#####################################################################################
# Calculate the mappable genome size by taking the full genome size and
# subtracting from it the length of all of the gaps. 
# For rhemac2.0 92% is non-gap.
# For rheMac8.0 only 88% is non-gap.
#####################################################################################    
WHICHBUILD = 'RHEMAC2'
if WHICHBUILD == 'RHEMAC2':
    chromsizefile = open('/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/RheMacChromSizes.txt')
    gapfile = open( '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/gap.txt')
elif WHICHBUILD == 'RHEMAC8':
    chromsizefile = open('/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF8/rheMac8.chrom.sizes.txt')
    gapfile = open( '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac8/rheMac8.gap.txt')
else:
    raise ValueError('Unknown reference genome build')
    
chrsizedict = {}
while True:
    line = chromsizefile.readline()
    linelist = line.split()
    if line == '':
        break
    if linelist[0][:4] == 'chrU':
        continue
    chrsizedict[linelist[0]] = int(linelist[1])
chromsizefile.close()
totalgenomesize = sum(chrsizedict.values())

# Now getting the gaps
gapdict={}
while True:
    line = gapfile.readline()
    linelist = line.split()
    if line == '':
        break
    if linelist[0] in gapdict.keys():
        gapdict[linelist[0]] += int(linelist[5])
    else:
        gapdict[linelist[0]] = int(linelist[5])
gapfile.close()
totalgapsize = sum([k for v,k in gapdict.items() if v in set(chrsizedict.keys())]) # Only looking at gaps in the standard chromosomes, avoiding "chrUn_blahblah".
mappablegenomesize = totalgenomesize-totalgapsize
mappablegenomesize/totalgenomesize

#%%
#####################################################################################
# Finding out which genes (an how many) are near gaps
#####################################################################################

import HTSeq
from matplotlib import pyplot
import numpy as np

gappos = GetGapPos()
gtffile = HTSeq.GFF_Reader("/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf")

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

genehitprctiles = [50, 60, 70, 80, 90]
cutoffvalues = np.percentile(RNAdict.values(),genehitprctiles)
cutoffvalues = np.append(cutoffvalues, maxreads)

# Getting the ranks of the RNA expression levels (using raw read counts)
RNAranks = {}
for a, b in RNAdict.iteritems():
    for j in range(len(cutoffvalues)):
       if a not in RNAranks and b <= cutoffvalues[j] :
           RNAranks[a] = j # ranks of gene expression
           break

gapdict = {}
for feature in gtffile:
    if feature.type == "exon" and feature.attr["exon_number"] == "1":
         p = feature.iv.start_d_as_pos # start_d_as_pos = The position of the TSS irrespective of strand
         if (p.start > halfwinwidth):
            window = HTSeq.GenomicInterval( p.chrom, p.pos - halfwinwidth, p.pos + halfwinwidth, "." )
            steps = [i for i in gappos[window].steps()]
            if len(steps) > 1:  # There should only be one "step" in c if there's no gap near the TSS
                gapdict[feature.name] = 1 # 1 mean "there is a gap"
            else:
                gapdict[feature.name] = 0
ngenes = len(gapdict.values())
ngeneswithgaps = sum(gapdict.values())
print("fraction of genes culled because of gaps near TSS: {}",format(ngeneswithgaps/float(ngenes)))  # Holy cow! 61% of genes have a gap in within 6 kb of the TSS
# Spot checking
genenames = list(gapdict.keys())

somegenenames = [genenames[int(i)] for i in np.rint(np.linspace(0,ngenes-1,10))]
somegapbools = [gapdict[k] for k in somegenenames]
print(zip(somegenenames,somegapbools))
print('\n'.join('{}: {}'.format(i,j) for i,j in zip(somegenenames,somegapbools)))