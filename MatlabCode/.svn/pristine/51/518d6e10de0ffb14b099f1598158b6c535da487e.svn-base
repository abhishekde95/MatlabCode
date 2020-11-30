
########################################################################
# Looking at the spatial relationship between TSSs and gaps.
# Plan is to get the mean distance between each TSS and the nearest
# gap. Then repeat this, scrambling the positions of the gaps.
# ######################################################################
import numpy as np
import HTSeq
import math
from matplotlib import pyplot
from itertools import compress
import copy

# Got from the web
# http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

# A useful representation for gappos might be a dictionary each key of which is a chromosome and each
# key is an ndarray of gap starts and stops? Think about how we're going to do the scrambling.
def GetGapPos(filename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/gap.txt'):
    '''Function to make a GenomicInterval object that knows where all the gaps in the rhemac2.0 sequence.
    Returns a GenomicInterval.'''
    import HTSeq
    gapfile = open(filename)
    gapdict = dict()
    while True:
        line = gapfile.readline()
        linelist = line.split()
        if line == '':
            break
        if linelist[0] in gapdict:
            newrow = np.fromiter(linelist[1:3], dtype=int)
            newrow.shape = (1, 2)
            gapdict[linelist[0]] = np.concatenate((gapdict[linelist[0]], newrow), axis=0)
        else:
            gapdict[linelist[0]] = np.expand_dims(np.fromiter(linelist[1:3], dtype = int), axis = 0)
    gapfile.close()
    for k in gapdict.keys():
        gapdict[k].sort(axis = 0) # gaps have to be sorted for find_nearest to work
    return(gapdict)

# Finds the shortest distance between a genomic location and a gap. Returns
# zero if the location is inside a gap.
# Input arguments:
#     tsslist: a list of lists. Each sublist consists of a chromosome (e.g. 'chr7') and position (e.g. 1151820)
#     gapdict: a dictionary. Keys are chromosomes (chr1-X, no MT). Values are nx2 ndarrays containing gap starts
# and stops.
def GetDistances(tsslist, gapdict):
    distances = np.zeros(shape=len(tsslist))
    counter = 0
    for i in tsslist:
        try:
            gaps = gapdict[i[0]]
            idx0 = np.searchsorted(gaps[:, 0], i[1], side="left") # nearest gap start
            idx1 = np.searchsorted(gaps[:, 1], i[1], side="left") # nearest gap stop
            #print([gaps[idx0,0], i[1], gaps[idx1,1]])
            if idx0 >= len(gaps):
                idx0 = idx0-1
            if idx1 >= len(gaps):
                idx1 = idx1-1
            if i[1] >= gaps[idx0,0] and i[1] <= gaps[idx1, 0]:
                distances[counter] = 0
                print('TSS inside gap: {}:{}'.format(i[0],i[1]))
            else:
                candidates = np.asarray([])
                if idx0 < len(gaps):
                    candidates = np.append(candidates, gaps[idx0,0])
                if idx0 < len(gaps):
                    candidates = np.append(candidates,gaps[idx1,1])
                if idx0 > 0:
                    candidates = np.append(candidates,gaps[idx0-1,0])
                if idx1 > 0:
                    candidates = np.append(candidates,gaps[idx1-1,1])
                d = i[1]-candidates
                distances[counter] = min(abs(d))
            counter += 1
        except KeyError:
            pass
        except IndexError:
            print(counter)
            print(idx0, idx1)
            print(candidates, d)
            print(len(tsslist), len(distances), len(gaps))
            break
    return distances

# Works but isn't that useful right now because some TSSs overlap with gaps
# naturally (!?) At least on chr1
def Overlap(tss, gaps):
    found_one = False
    for i in tss:
        idx = np.searchsorted(gaps[:,0], i, side="left")
        if gaps[idx,0] == i:
            found_one = True
        elif i > gaps[idx,0]: # (TSS is downstream of nearest gaps start)
            if i < gaps[idx,1]:
                found_one = True
        else: # i < gaps[idx,0] (TSS is upstream of nearest gap start)
            if idx > 0:
                if i >= gaps[idx-1,0] and i <= gaps[idx-1,1]:
                    found_one = True
        if found_one:
             break
    if found_one:
        return(True)
    else:
        return(False)

# This is going to be tricky. I want to move TSSs around within chromosomes preserving their interval
# distribution and making sure that none of them fall on a gap.
def ScrambleTssdict(tssdict):
    import random
    for chr in tssdict.keys():
         if chr in tssdict:
             tmptsslist = tssdict[chr]
             intervals = np.diff(np.asarray(tmptsslist))
             random.shuffle(intervals)
             tssdict[chr] = np.cumsum(np.insert(intervals, 0, tmptsslist[0]))
    return tssdict

         #counter = 0
         #while(Overlap(np.cumsum(np.insert(intervals,0,tmptsslist[0])),gapdict[chr])):
         #    random.shuffle(intervals)
         #    counter += 1
         #    print(counter)
# Still working on this ^

tsslist = list()
tssdict = dict()
gtffile = HTSeq.GFF_Reader("/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf")
for feature in gtffile:
    if feature.type == "exon" and feature.attr["exon_number"] == "1" and feature.iv.chrom[0:3] == 'chr' and not feature.iv.chrom[0:4] == 'chrM':
       p = feature.iv.start_d_as_pos  # start_d_as_pos = The position of the TSS irrespective of strand
       tsslist.append([p.chrom, p.pos]) # chrom, tsspos
       if p.chrom in tssdict:
           tssdict[p.chrom].append(p.pos)
       else:
           tssdict[p.chrom] = [p.pos,]
for k in tssdict:
    tssdict[k].sort()
preservedtssdict = copy.deepcopy(tssdict) # THIS ISN'T WORKING!

gapdict = GetGapPos()

# Just confirming with the genome browser that the TSS that I think is farthest from a gap
# really is that distance from a gap. Appears to check out.
#t = distances==np.max(distances)
#i = tsslist[list(compress(xrange(len(t)), t))[0]]
#find_nearest(np.ravel(gapdict[i[0]]), i[2])

# Now trying a permutation test to see if TSSs really are preferentially near gaps in rhemac2.

niter = 1000
data = np.zeros(shape=(niter,1))
for i in range(niter):
    if i == 0:
        #tmptsslist = tsslist
        tmptssdict = preservedtssdict
    else:
        tmptssdict = ScrambleTssdict(tssdict) # destructive modification
    tmptsslist = list()
    for k in tmptssdict.keys():
        for v in tmptssdict[k]:
            tmptsslist.append([k,v])
    #distances = np.zeros(shape=len(tsslist))
    distances = GetDistances(tmptsslist, gapdict)
    ntssinsidegaps = np.sum(distances == 0)
    #pyplot.plot(distances)
    mndist = np.mean(abs(distances))
    print("Iter {}: mean distances is {}. {} TSSs inside gaps.".format(i, mndist, ntssinsidegaps))
    data[i] = mndist

pyplot.hist(data[1:])
pyplot.plot(data[0],0,'m*')
pyplot.title('Distance of TSS to nearest gap in rhemac2')
pyplot.savefig('rhemac2_TSSvsgapMC.pdf')