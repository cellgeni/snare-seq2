#!/usr/bin/env python
import gzip
import argparse
import json
import sys

parser = argparse.ArgumentParser(description='Extracts and rearrange barcodes and umis into artificiall fastq and prints it in stdout.')
parser.add_argument("infile", type=str,default=None,help='input fastq.gz file')
#parser.add_argument("output", type=str,default=None,help='output fasta.gz file')
parser.add_argument("coord", type=str,default=None,help='coordinates of pieces to combine in the following form: "0-10,22-30". Zero-based and left inclusive and right non-inclusive')
parser.add_argument("whitelists", type=str,default=None,help='list of paths (separated by comma) to whitelist files, one per interval.')

args = parser.parse_args()

#args = parser.parse_args(['data/GBM_spatial_m13050790/45611_1_5_R2_001.top3.fastq.gz','32-40,70-78,22-32','actions/DBiT-seq/Datasets/barcode8.txt,actions/DBiT-seq/Datasets/barcode8.txt,'])


def oneMismHash(ss,a=['A','T','G','C']):
  r = {w:w for w in ss}
  for s in ss:
    for i in range(0,len(s)):
      for l in a:
        if l != s[i]:
          m = s[0:i] + l + s[(i+1):]
          r[m] = s
  return r


# read whitelists
whitelists0 = []
whitelists1 = []
stat = {'general':{'nreads':0,'wl0':0,'wl1':0},'perbarcode':[]}
for fn in args.whitelists.split(','):
  stat['perbarcode'].append({'wl0':0,'wl1':0})
  if fn != '':
    with open(fn) as f:
      bcs = f.read().splitlines()
      whitelists0.append(set(bcs))
      whitelists1.append(oneMismHash(bcs))
  else:
    whitelists0.append({})
    whitelists1.append({})
    
# prepare coordinates
coordinates = [[int(i) for i in s.split('-')] for s in args.coord.split(',')]

# process
fin = gzip.open(args.infile,'rt')
fout = sys.stdout#open(args.output,'wt')#gzip.open(args.output,'wt')

while True:
  read = [fin.readline() for i in range(0,4)]
  if read[0] == '':
    break
  stat['general']['nreads'] += 1
  # fix qual
  q = ''
  for c in coordinates:
    q = q + read[3][c[0]:c[1]]
  read[3] = q+"\n"
  # fix read
  r = ''
  good0 = True
  good1 = True
  for i in range(0,len(coordinates)):
    c = coordinates[i]
    f = read[1][c[0]:c[1]]
    if len(whitelists1[i]) > 0:
      if f in whitelists0[i]:
        stat['perbarcode'][i]['wl0'] += 1
        stat['perbarcode'][i]['wl1'] += 1
      else:
        good0 = False
        if f in whitelists1[i]:
          f = whitelists1[i][f]
          stat['perbarcode'][i]['wl1'] += 1
        else:
          good1 = False
    else:
      stat['perbarcode'][i]['wl0'] += 1
      stat['perbarcode'][i]['wl1'] += 1
    r = r + f
  if good0:
    stat['general']['wl0'] += 1
  if good1:
    stat['general']['wl1'] += 1
  read[1] = r+"\n"
  fout.writelines(read)

fin.close()
#fout.close()

json.dump(stat, sys.stderr)

