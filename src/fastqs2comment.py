#!/usr/bin/env python
import gzip
import sys
import argparse

parser = argparse.ArgumentParser(description='sets sequence from one fastq as comment in another fastq and prints resulting fq into std.')
parser.add_argument("fastq", type=str,default=None,help='input fastq.gz file')
parser.add_argument("comment", type=str,default=None,help='input fastq.gz file with sequence to be used as comment')
parser.add_argument("prefix", type=str,default=None,help='prefix to add to sequence (use "CB:Z:" to add cell barcode)')

args = parser.parse_args()

fin = gzip.open(args.fastq,'rt')
cin = gzip.open(args.comment,'rt')


while True:
  read = [fin.readline() for i in range(0,4)]
  if read[0] == '':
    break
  comment = [cin.readline() for i in range(0,4)]
  read[0] = read[0].split(' ')[0].strip() + ' ' + args.prefix + comment[1]
  sys.stdout.writelines(read)


fin.close()
cin.close()

