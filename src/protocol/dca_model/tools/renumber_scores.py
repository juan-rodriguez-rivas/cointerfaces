#!/usr/bin/env python
import sys 
import os
import numpy as np

def get_command():
    import argparse
    parser = argparse.ArgumentParser(description='Renumber scores to their HMM numbering and apply APC')
    parser.add_argument('-d','--fdat',type=str, help="dat file")
    parser.add_argument('-s','--fscores',type=str, help="scores file")
    parser.add_argument('-o','--out_root',type=str, help="output root",default='./')
    args = parser.parse_args()

    fdat =  args.fdat
    fscores =  args.fscores
    out_root = args.out_root

    if not fdat or not fscores: 
        parser.print_help()
        exit()

    return(fdat,fscores,out_root)

# read args and options from command line
fdat,fscores,out_root = get_command()

# read from dat
nlines=0
for line in open(fdat,'r'): 
    nlines+=1
    nv1,nv2,ns=map(int, line.split()[:3]) # read the first 3 fields only 
    break # read the first line only

# retrieve label
label='.'.join(os.path.basename(fscores).split(".")[:-1])

# read scores 
with open(fscores) as f: lines=[line for line in f]
nlines=sum([1 for line in lines if len(line.split())==3])
if int(np.rint(0.5*(np.sqrt(8.0*nlines)+1.0))) != nv1+nv2: 
    sys.stderr.write("ERROR: check num. of positions")
    exit()
mat=np.zeros((nv1+nv2,nv1+nv2))
for line in lines: 
    a,b,s=line.split()
    a=int(a)
    b=int(b)
    mat[a-1,b-1]=float(s)
    mat[b-1,a-1]=float(s)

# compute averages over rows/columns

sm1=np.mean(mat[:nv1,:nv1])
cm1=np.mean(mat[:nv1,:nv1],axis=0) 
sm2=np.mean(mat[nv1:,nv1:])
cm2=np.mean(mat[nv1:,nv1:],axis=0)
sm12=np.mean(mat[:nv1,nv1:])
cm12=np.mean(mat[:nv1,nv1:],axis=0)
rm12=np.mean(mat[:nv1,nv1:],axis=1)

# correct for APC a la Baker
mat_apc=np.copy(mat)
for a in range(nv1): 
    for b in range(nv1): 
        apc = cm1[a]*cm1[b]/sm1
        mat_apc[a,b]=mat[a,b]-apc
for a in range(nv2): 
    for b in range(nv2): 
        apc = cm2[a]*cm2[b]/sm2
        mat_apc[nv1+a,nv1+b]=mat[nv1+a,nv1+b]-apc
for a in range(nv1): 
    for b in range(nv2): 
        apc = rm12[a]*cm12[b]/sm12
        mat_apc[a,nv1+b]=mat[a,nv1+b]-apc
        mat_apc[nv1+b,a]=mat[nv1+b,a]-apc

u1=open(out_root+'.scores.intra1','w')
u2=open(out_root+'.scores.intra2','w')
uint=open(out_root+'.scores.inter','w')

for a in range(nv1+nv2-1): 
    for b in range(a+1,nv1+nv2):
        if b <= a: continue
        inda=a+1
        indb=b+1
        if inda <= nv1 and indb <= nv1: 
            # intra-domain1
            sp=mat_apc[a,b]
            u1.write(' '.join([str(x) for x in [inda,indb,sp]])+'\n')
        elif inda <= nv1 and indb > nv1: 
            # inter-domain1-domain2
            sp=mat_apc[a,b]
            uint.write(' '.join([str(x) for x in [inda,indb-nv1,sp]])+'\n')
        elif inda > nv1 and indb > nv1: 
            # intra-domain2
            sp=mat_apc[a,b]
            u2.write(' '.join([str(x) for x in [inda-nv1,indb-nv1,sp]])+'\n')

