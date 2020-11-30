#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:05:57 2017

@author: horwitzlab
"""

# This *has* to be run as iPython, not Python because of the "magic" commands (prefixed with "!")

# Remember,  when saving BED files from Excel save as "windows formatted text"
# Also important: each line should end with a linefeed ("\n") not a carriage return 
# ("\r") nor a combination (e.g. "\n\r"). This can be done in TextEdit with a search and replace
# (click on the magnifying glass)
from matplotlib import pyplot
import numpy as np
import itertools
import os
import dill
import subprocess
import MapReduceHackery as MRH
import pybedtools


macs2out_filename = MRH.findfile('Z12131_all_peaks_summits.bed')
#macs2out_filename = MRH.findfile('A13062_V1_nuc_peaks_summits.bed')
blacklistfilename = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2/rhemac2blacklist.bed'
blackregions = pybedtools.BedTool(blacklistfilename)
blackregions.head()


# Trying to call annotatePeaks.pl from here.
# First need to make a BED file that HOME can deal with, with columns:
# 1) unique peak identifier
# 2) chromosome
# 3) starting position
# 4) ending position
# 5) strand (using "1" in all cases since ATAC-seq is unstranded

if macs2out_filename is None:
    my_error = IOError("Cannot find mac2 output file")
    raise my_error
macs2out = pybedtools.BedTool(macs2out_filename)

# Making a BEDfile that annotate.pl likes and saving it with the name "tmpfileforhomer"
data = []
counter = 0
for i in macs2out:
    data.append((str(counter),str(i[0]),str(i[1]),str(i[2]),'1'))
    counter +=1
tmpfileforhomer = open('tmpfileforhomer','w')
for i in data:
    for j in i:
        tmpfileforhomer.write(j)
        tmpfileforhomer.write("\t")
    tmpfileforhomer.write("\n")
tmpfileforhomer.close()
var = 'tmpfileforhomer /Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/FASTA_ETC/rheMac2 -gtf /Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf -gid -organism Mmulatta'
annotatePeaks_call_string = 'annotatePeaks.pl'

annotatePeaks_dir ="/Users/horwitzlab/Desktop/SEQ_ANALYSIS/homer/bin"
%env PATH=/Users/horwitzlab/anaconda/bin:/usr/bin:/bin:/usr/sbin:/sbin:{annotatePeaks_dir}
# Below, running annotatePeaks.pl from within IPython!
# Not sure how I'm going to retrieve the data
tempfilename = 'tmpfile.tmp'
!{annotatePeaks_call_string+' '+var+' > '+ tempfilename}

# Reading in an annotation file from the Homer program, annotatePeaks.pl
annotationfilename = tempfilename
annotationfile = open(annotationfilename)
data = []
while True:
    line = annotationfile.readline()
    if line == '':
       break
    if not line[0].isalpha():
       linelist = line.split()
       if not linelist[-2].isalpha():
          data.append((linelist[1],linelist[2],linelist[3],linelist[-2],linelist[-1]))
annotationfile.close()
os.remove(tempfilename)
# Getting rid of blacklisted peaks
annotated_peaks = pybedtools.BedTool(data)
a_and_b = annotated_peaks.intersect(blackregions)
real_peaks = annotated_peaks-a_and_b

# Getting distances between peaks and nearest TSSs
# Note: According to the Homer documentation "distance to TSS" is signed such that
# distances < 0 are upstream and distances > 0 are downstream.
data = []
for i in real_peaks:
    data.append(int(i[3]))

data = [i if v else 0.1 for (i,v) in zip(data, np.abs(data) != 0)] # replacing "0" with ".1" to avoid log(0) errors (pythonically)
logdists = np.transpose(np.log10(np.abs(np.asmatrix(data))))
logdists = np.multiply(np.transpose(np.asmatrix(np.sign(data))),logdists)
pyplot.figure()
pyplot.hist(logdists,bins = 200)
pyplot.title(macs2out_filename[macs2out_filename.rfind(os.path.sep)+1:])
pyplot.xlabel('log10 distance from TSS')
pyplot.ylabel('count')
pyplot.xticks(np.arange(-7,8,1))
#pyplot.yscale('log')
#%%
# Getting the names of the genes that are near significant peaks
# (within particular distances from genes)

gtfname = '/Users/horwitzlab/Desktop/SEQ_ANALYSIS/genomes/GTF/MmulattaGTF/Macaca_mulatta.MMUL_1.83.chr.gtf'
genes = pybedtools.BedTool(gtfname)

ensembldict = dict()
reverseensembldict = dict()
# Making a dictionary that will map ENSEMBL codes onto gene names
# takes a surpringly long time (about a minute)
for i in genes:
    if i.attrs.has_key('gene_id') and i.attrs.has_key('gene_name'):
        if not ensembldict.has_key(i.attrs['gene_id']):
            ensembldict[i.attrs['gene_id']] = i.attrs['gene_name']
        if not reverseensembldict.has_key(i.attrs['gene_name']):
            reverseensembldict[i.attrs['gene_name']] = i.attrs['gene_id']
#%%  
query = ''
while (query != 'exit'):    
    query = raw_input('Enter a gene name (type "exit" to quit): ')
    if reverseensembldict.has_key(query):
       print ('ENSEMBL code is: {}'.format(reverseensembldict[query]))
       for i in real_peaks:
           if str(i.fields[-1]) == reverseensembldict[query]:
               print(i)
    else:
        if query != 'exit':
            print('Could not find gene named {}\n'.format(query))
#%%
# Loading in some RNA seq for an automated way of finding genes that are highly expressed in the macaque brain
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
cutoffvalue = np.percentile(RNAdict.values(),90)
whichgenes = set([i for i in RNAdict.keys() if RNAdict[i] > cutoffvalue])
for i in real_peaks:
    ensemblid = str(i.fields[-1])
    if ensemblid in whichgenes:
        if ensembldict.has_key(ensemblid):
            print('{} {}'.format(ensembldict[ensemblid],i))
        else:
            print(i)