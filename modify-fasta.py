#!/usr/bin/env python3.8
# -*- coding: utf-8 -*-

"""
@author: Dilek Koptekin
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from pybedtools import BedTool

def get_all(dfBed, chrom, chrLen):
    pstart = 0
    posList = []
    for index, row in dfBed.iterrows():
        pend = int(row['posEnd'])
        posList.append([chrom, pstart, pend-1, "_".join(map(str, [chrom, pstart, pend-1]))])
        posList.append([chrom, pend-1,pend, "_".join(map(str, [chrom, pend-1, pend]))])
        pstart = pend
    if int(pstart) < chrLen:
        pend = chrLen
        posList.append([chrom, pstart, pend, "_".join(map(str, [chrom, pstart, pend]))])
        df = pd.DataFrame(posList, columns=['chr', 'posStart', 'posEnd', 'snpID'])
    return df[df['posStart'] != df['posEnd']]

def convert2fasta(allSeq, snpFile, nFa):
    allSeq.loc[allSeq.snpID.isin(snpFile.snpID), ['seq']] = snpFile[['allele2']].values
    fa2 = "".join(allSeq['seq'])
    fa2_parts = [fa2[i:i + 60] for i in range(0, len(fa2), 60)]
    if nFa == 1:
        return fa2_parts
    if nFa == 2:
        allSeq.loc[allSeq.snpID.isin(snpFile.snpID), ['seq']] = snpFile[['allele1']].values
        fa1 = "".join(allSeq['seq'])
        fa1_parts = [fa1[i:i + 60] for i in range(0, len(fa1), 60)]
        return fa1_parts, fa2_parts

        sys.exit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='modify-fasta')
    parser.add_argument('-fi','--fasta', type=str,  dest='fasta_file', required=True, help='')
    parser.add_argument('-bed', type=str,  dest='bed_file', required=True, help='')
    parser.add_argument('-o', '--output', dest='output_file', help='')
    parser.add_argument('-n', type=int, dest='n_fa', default=1)
    args = parser.parse_args()

    fi = args.fasta_file
    bed = args.bed_file
    n_fa = args.n_fa

    if args.output_file is not None:
        output = args.output_file
    else:
        output = os.path.splitext(os.path.basename(bed))[0]

    fai_path = os.path.join(fi + ".fai")

    if os.path.exists(fai_path):
        fai = pd.read_csv(fai_path, delim_whitespace=True, header=None,usecols=[0,1], names = ["chr", "length"])
        fai = fai.set_index('chr')['length'].to_dict()
    else:
        print(f'{fai_path} not found \n'
              f'Please index the fasta file by using "samtools faidx"')

    chrs = np.unique(np.loadtxt(bed, dtype=str, usecols=0))
    for i in chrs:
        if fai.get(i) == None:
            print('The headers in the input FASTA file do not match the chromosome column in the BED file.')
            sys.exit()
    n_chr = len(chrs)

    if n_fa == 1:
        fasta2 = open(os.path.join(output + '.fa'), 'w')
    elif n_fa == 2:
        fasta1 = open(os.path.join(output + '.1.fa'), 'w')
        fasta2 = open(os.path.join(output + '.2.fa'), 'w')
    else:
        print("Number of output fasta could be 1 or 2")
        sys.exit()

    snp = pd.read_csv(bed, delim_whitespace=True, names=['chr', 'posEnd', 'allele1', 'allele2'],
                      dtype={'chr':'str', 'posEnd':'int', 'allele1':'str', 'allele2':'str'})

    snp.insert(1, 'posStart', snp['posEnd']-1)
    snp['snpID'] = snp['chr'] + "_" + snp['posStart'].astype(str) + "_" + snp['posEnd'].astype(str)

    for chrom in chrs:
        chr_length = fai.get(chrom)
        snp_chr = snp[snp.chr == chrom]
        df = get_all(snp_chr, chrom, chr_length)
        if df[df['posStart'] > df['posEnd']].size > 0:
            print('The input bed file must be sorted by chr and position and include only biallelic SNPs.')
            sys.exit()
        all_seq_chr = pd.read_csv(BedTool.from_dataframe(df).sequence(fi=fi, tab=True, name=True).seqfn, sep="\t", names=["snpID","seq"])
        if n_fa == 1:
            fa2 = convert2fasta(all_seq_chr, snp_chr, n_fa)
            fasta2.writelines(f">{chrom}\n")
            fasta2.writelines("\n".join(fa2))
            fasta2.writelines("\n")
        if n_fa == 2:
            fa1, fa2 = convert2fasta(all_seq_chr, snp_chr, n_fa)
            fasta1.writelines(f">{chrom}\n")
            fasta1.writelines("\n".join(fa1))
            fasta1.writelines("\n")
            fasta2.writelines(f">{chrom}\n")
            fasta2.writelines("\n".join(fa2))
            fasta2.writelines("\n")
    if n_fa == 1:
        fasta2.close()
    elif n_fa == 2:
        fasta1.close()
        fasta2.close()
    print("Finished")
