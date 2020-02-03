#!/usr/bin/env python
import sys 
import os
import numpy as np
from Bio import SeqIO

# list of natural aminoacids
alphabet = [
'-',
'A','C','D','E',
'F','G','H','I',
'K','L','M','N',
'P','Q','R','S',
'T','V','W','Y',
]

# add 1 to the indices (gaps will have index 1)
AAMap = {alphabet[i]: i+1 for i in range(len(alphabet))}

def get_command():
    import argparse
    parser = argparse.ArgumentParser(description='clustering sequences')
    parser.add_argument('-1','--fasta1',type=str, help="first domain - FASTA file")
    parser.add_argument('-2','--fasta2',type=str, help="second domain - FASTA file")
    parser.add_argument('-o','--out_file',type=str, help="output file",default='./')
    parser.add_argument('-r','--n_shuf',type=int, help="n. of random inputs",default=0)
    args = parser.parse_args()

    ffasta1 =  args.fasta1
    ffasta2 =  args.fasta2
    out_file = args.out_file
    n_shuf = args.n_shuf

    if not ffasta1 or not ffasta2: 
        parser.print_help()
        exit()

    return(ffasta1,ffasta2,out_file,n_shuf)

# read args and options from command line
ffasta1,ffasta2,out_file,n_shuf  = get_command()

# parse the MSA files
try: 
    fasta_sequences_1 = list(SeqIO.parse(open(ffasta1),'fasta'))
except:
    sys.stderr.write("check FASTA file: "+ffasta1+"\n")
    exit()

try: 
    fasta_sequences_2 = list(SeqIO.parse(open(ffasta2),'fasta'))
except:
    sys.stderr.write("check FASTA file: "+ffasta2+"\n")
    exit()

if len(fasta_sequences_1) != len(fasta_sequences_2): 
    sys.stderr.write("n. of sequences in the MSA differ\n")
    exit()
else: 
    ns0 = len(fasta_sequences_1)
    
# msa is a list of sequences
msa1=[]
msa2=[]
msa=[]
# headers is a list of header
headers1=[]
headers2=[]

# loop over the joint alignment
for i in range(ns0):
    name1, sequence1 = fasta_sequences_1[i].id, fasta_sequences_1[i].seq.tostring()
    name2, sequence2 = fasta_sequences_2[i].id, fasta_sequences_2[i].seq.tostring()
    if i == 0: 
        n1 = len(sequence1)
        n2 = len(sequence2)
        ntot = n1+n2
    # filter out unnatural aminoacids
    if np.sum([1 for aa in sequence1 if not aa in alphabet]) > 0: continue
    if np.sum([1 for aa in sequence2 if not aa in alphabet]) > 0: continue
    headers1.append(name1)
    headers2.append(name2)
    msa1.append([ x for x in sequence1])
    msa2.append([ x for x in sequence2])
    msa.append([ x for x in sequence1+sequence2])

# switch to a numpy matrix for slicing
mat1=np.asarray(msa1)
mat2=np.asarray(msa2)
mat=np.asarray(msa)

# n. of sequences
ns = len(mat[:,0])
# n. of positions in domain1
nv1 = len(mat1[0,:])
# n. of positions in domain2
nv2 = len(mat2[0,:])
# total n. of positions 
nv = len(mat[0,:])

# dump a temp file
fprm=open(out_file + '.dim','w')
fprm.write("%i %i %i #np1 np2 nseq \n" % (nv1, nv2, ns))


if n_shuf == 0: 
    with open(out_file,'w') as f: 
        # dump input file(s) for mpl
        for s in range(ns): 
            f.write(' '.join('%2i' % AAMap[x] for x in mat[s,:])+'\n')
else:
    for i in range(n_shuf): 
        # shuffled indices
        inds=np.arange(ns)
        np.random.shuffle(inds)
        mat[:,nv1:]=mat2[inds,:]
        with open(out_file+'.shuf_'+str(i),'w') as f: 
            # dump input file(s) for mpl
            for s in range(ns): 
                f.write(' '.join('%2i' % AAMap[x] for x in mat[s,:])+'\n')
        
