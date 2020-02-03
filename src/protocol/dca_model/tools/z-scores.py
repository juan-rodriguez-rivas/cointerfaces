#!/usr/bin/env python 
# Copyright (c) 2015 Simone Marsili

import sys
import numpy as np

description = """
=====================================
z-scores - convert scores to z-scores 
=====================================
"""

def gumbel_cdf(x,n):
    if n > 1:
        logn = np.log(n)
        slogn = np.sqrt(2.0 * logn)
        mu = slogn - 0.5 * (np.log(logn) + np.log(4.0 * np.pi)) / slogn
        beta = 1.0 / slogn
        z = (x - mu) / beta
        return np.exp(-np.exp(-z))
    else:
        return 0.0

def get_command(description):
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=description, epilog=" ")
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument('-f',"--field", type=int, help="target field", required=True)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'), help='Input file name', default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), help='Output file name', default=sys.stdout)
    parser.add_argument("--mode", type=str, help="standardization option (std, mad)", default='mad')
    parser.add_argument("--pvalue", action='store_true', help='compute the p-value for i-th largest z-score to be the largest among n - i + 1 samples from a standard normal', default=False)
    parser.add_argument("--sort", action='store_true', help='sort z-scores for printing', default=False)
    parser.add_argument("--print_top", type=str, help='dump entry if z-score > 4 or p-value < 0.05 ', default=None)

    args = parser.parse_args()
    fin = args.infile
    fout = args.outfile
    field = args.field
    mode = args.mode
    pvalue = args.pvalue
    sort = args.sort
    print_top = args.print_top

    # check field value
    line = fin.readline()
    nfields = len(line.split())
    fin.seek(0) # rewind
    if field == 0: 
        sys.stderr.write("0 is not a valid field\n")
        exit()
    if field > nfields: 
        sys.stderr.write("target field ({!s}) is larger than the total number of fields ({!s})\n".format(field, nfields))
        exit()
    
    # check mode value
    if not mode in ['mad','std']: 
        sys.stderr.write("unkown standardization option ({!s})\n".format(mode))
        exit()

    # check print_top
    if print_top and not print_top in ['z-score', 'p-value']: 
        sys.stderr.write('possible arguments for --print_top are "z-score" and "p-value"\n')
        exit()
    if print_top=='p-value' and not pvalue: 
        # switch pvalue to True
        pvalue = True

    return(fin, fout, field, mode, pvalue, sort, print_top)

fin, fout, field, mode, pvalue, sort, print_top = get_command(description)

lines = [line for line in fin if line[0]!="#"]
nlines = len(lines)
try: 
    scores = np.asarray([float(line.split()[field-1]) for line in lines])
except:
    sys.stderr.write("check field {:d} for non-numerical values\n".format(field))
    exit()

if mode == 'mad': 
    center = np.median(scores)
    scores = scores - center
    sigma = 1.4826 * np.median(np.abs(scores))
    scores = scores / sigma 
elif mode == 'std': 
    center = np.mean(scores)
    scores = scores - center
    sigma = np.std(scores)
    scores = scores / sigma 
    
sorted_indices = np.argsort(scores)[::-1]
ranks = {sorted_indices[i]: i for i in range(nlines)}

for l in range(nlines):
    if sort:
        ind = sorted_indices[l]
        rank = l
    else:
        ind = l
        rank = ranks[l]
    fi = lines[ind].split()
    s = scores[ind]
    fi[field-1] = "{:.5f}".format(s)
    if pvalue:
        pval = 1.0 - gumbel_cdf(s,nlines-rank)
        fi.append("{:.2e}".format(pval))
    if print_top == 'p-value': 
        if pval >= 0.05: continue
    elif print_top == 'z-score': 
        if s < 4: continue
    print('\t'.join(fi))
    
