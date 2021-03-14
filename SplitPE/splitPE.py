#!/usr/bin/env python3

##############################################################################
# FUNCS

bases=[x for x in "ACGT"]
def getMM1Barcodes(barcodes):
    barcodesm1=copy.deepcopy(barcodes)
    for bcstr in barcodes.keys():
        bcs=[x for x in bcstr]
        for i in range(len(bcs)):
            bcst=bcs[:]
            for bi in bases:
                if bcs[i]!=bi:
                    bcst[i] = bi
                    barcodesm1["".join(bcst)]=barcodes[bcstr]
    return barcodesm1

def writeFastq(fp,rr,bcid):
    print(rr[0],rr[3],bcid,file=fp)
    print(rr[1],file=fp)
    print("+",file=fp)
    print(rr[2],file=fp)

#
##############################################################################

import sys
import copy
import os
import gzip

try:
    import pyfastx
except ModuleNotFoundError:
    print("\n   Need to have pyfastx installed\n")
    print("\n   https://github.com/lmdu/pyfastx\n")
    sys.exit(1)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-b",help="Run in block mode [-b m:n block m of n]")
parser.add_argument("Barcodes",help="Barcode File")
parser.add_argument("R1",help="R1 Fastq file")
parser.add_argument("R2",help="R2 Fastq file",nargs="?",default="")
args=parser.parse_args()

if args.R2=="":
    args.R2=args.R1.replace("_R1_","_R2_")

barcodes=list()
with open(args.Barcodes) as fp:
    for line in fp:
        (bc,seq)=line.strip().split("\t")
        barcodes.append((seq,bc))
barcodes=dict(barcodes)

#
# Do a subblock
#

if args.b:
    subBarcodes=dict()
    (m,n)=[int(x) for x in args.b.split(":")]
    print(m,n)
    m=m-1
    numBarcodes=len(barcodes.values())
    for bi in barcodes.keys():
        bid=int(barcodes[bi][1:])
        if(bid % n == m):
            subBarcodes[bi]=barcodes[bi]
    barcodes=subBarcodes

barcodes1MM=getMM1Barcodes(barcodes)

base=os.path.basename(args.R1).split("_R1_")[0]
obase=os.path.join("splits",base)
outFiles=dict()
for bcid in barcodes.values():
    odir=os.path.join(obase,bcid)
    os.makedirs(odir, exist_ok=True)
    splitR1File=os.path.join(odir,base+"__"+bcid+"__R1_000.fastq.gz")
    splitR2File=os.path.join(odir,base+"__"+bcid+"__R2_000.fastq.gz")
    fpR1=gzip.open(splitR1File,"wt")
    fpR2=gzip.open(splitR2File,"wt")
    outFiles[bcid]=(fpR1,fpR2)


import pyfastx

fq1=pyfastx.Fastx(args.R1)
fq2=pyfastx.Fastx(args.R2)

BCSTART=6
BCLEN=len("TTGTCAAGCAG")
BCEND=BCSTART+BCLEN


import time
start = time.time()

i=0
LIMIT=277684840
for r1,r2 in zip(fq1,fq2):
    i=i+1
    if i % 1000000 == 0:
        delta=time.time()-start
        minutes=int(delta/60)
        print(100.0*i/LIMIT, str(minutes)+":"+"%.3f" % (delta-minutes*60), i)
    # if i % LIMIT == 0:
    #     sys.exit()
    bc1=r1[1][BCSTART:BCEND]
    bc2=r2[1][BCSTART:BCEND]
    bcid1=barcodes1MM.get(bc1,"UNK")
    bcid2=barcodes1MM.get(bc2,"UNK")
    if bcid1!="UNK" and bcid1==bcid2:
        writeFastq(outFiles[bcid1][0],r1,bcid1+":"+bc1)
        writeFastq(outFiles[bcid1][1],r2,bcid2+":"+bc2)
