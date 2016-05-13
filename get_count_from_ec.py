#!/usr/bin/python
# usage: python get_count_from_ec.py normal_count.txt mgm_function_seed ec > output
import sys
count = {}
count_read = open(sys.argv[1],'r')
for n,line in enumerate(count_read):
    if n == 0:
        print line,
        continue
    spl = line.strip().split('\t')
    count[spl[0]] = line.strip()

anno_read = open(sys.argv[2],'r')
ec = sys.argv[3]
for n,line in enumerate(anno_read):
    #skip the first row
    if n == 0:
        continue
    spl = line.strip().split('\t')
    #skip the last row
    if (len(spl)<12):
        continue
    if (ec in spl[12]):
        #print line
        #get count
        #get id
        ids = '_'.join(spl[0].split('|')[1].split('_')[:-3])
        print count[ids]
