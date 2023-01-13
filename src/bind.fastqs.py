#!/usr/bin/env python
import gzip
import sys

infiles=[]
for fn in sys.argv[1:]:
  infiles.append(gzip.open(fn,'rt'))


while True:
  read = [infiles[0].readline() for i in range(0,4)]
  if read[0] == '':
    break
  for inf in infiles[1:]:
    read1 = [inf.readline() for i in range(0,4)]
    read[1] = read[1].rstrip()+read1[1]
    read[3] = read[3].rstrip()+read1[3]
  sys.stdout.writelines(read)

for inf in infiles:
  inf.close()

