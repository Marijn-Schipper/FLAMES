import pandas as pd
import os
import subprocess
import time
from sys import argv
import argparse
import gzip
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--filedir', required=True, help="Path to input directory where the file GenomicRiskLoci.txt is located.")
    parser.add_argument('--sample_size', required=True, help="Sample size")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    filedir = args.filedir
    sample_size = args.sample_size

    loci = pd.read_csv(os.path.join(filedir, 'GenomicRiskLoci.txt', sep='\t')) #make sure that the path is set correctly

    if not os.path.exists("loci"):
        os.makedirs("loci")
        
    # make the header.txt file
    # check if the first line is a header that begins with #
    with gzip.open(os.path.join(filedir, "input.snps.gz"), 'rt') as f:
        first_line = f.readline()
        if not first_line.startswith("#"):
            print("input.snps.gz does not have a valid header line starting with #.")
            sys.exit(1)
        else:
            header = first_line[1:].rstrip("\n").split("\t")
            print("\t".join(header), file=open(os.path.join(filedir, "header.txt"), 'w'))
    
    for index, row in loci.iterrows():
        locname = row['GenomicLocus']
        
        chrom = row["chr"]
        start = row["start"]
        end = row["end"]
        locname = f'{chrom}:{start}-{end}'
        if os.path.exists(f"locus_{locname}.susie.finemapped"):
            print(f"locus_{locname}.susie.finemapped already exists, skipping")
            continue

        cmd = f'cat header.txt > loci/locus_{locname}.txt' #from fuma: chr    bp      A2       A1   rsID    p       beta
        process = subprocess.Popen([cmd], close_fds=True, shell=True)
        process.wait()

        if start == end:
            start = start - 25000
            end = end + 25000

        cmd = f"tabix input.snps.gz {chrom}:{start}-{end} >> loci/locus_{locname}.txt"
        process = subprocess.Popen([cmd], close_fds=True, shell=True)
        process.wait()
        time.sleep(2.5)
        
    
        try:
            cmd = f'python /polyfun/munge_polyfun_sumstats.py --sumstats loci/locus_{locname}.txt --out loci/locus_{locname}.txt.pq --n {sample_size}'
            print(cmd)
            process = subprocess.Popen([cmd], close_fds=True, shell=True)
            process.wait()
            time.sleep(2.5)

            cmd = f'python /polyfun/finemapper.py --sumstats loci/locus_{locname}.txt.pq --method susie --n {sample_size} --out loci/locus_{locname}.susie --max-num-causal 1 --chr {chrom} --start {start} --end {end} --non-funct'
            process = subprocess.Popen([cmd], close_fds=True, shell=True)
            print(cmd)
            process.wait()
            time.sleep(2.5)
            
            finemapped = pd.read_csv(f"loci/locus_{locname}.susie", sep="\t")
            print(finemapped)
            print(max(finemapped['PIP']))
            finemapped = finemapped.sort_values("PIP", ascending=True)
            if len(finemapped) > 1:
                finemapped['cumsum'] = finemapped['PIP'].cumsum()
                finemapped = finemapped[finemapped['cumsum'] > 0.05]
            finemapped['SNP'] = finemapped.apply(lambda row: f"{row['CHR']}:{row['BP']}:{row['A1']}:{row['A2']}", axis=1)
            finemapped = finemapped[['SNP', 'PIP']]
            finemapped = finemapped.sort_values("PIP", ascending=False)
            if len(finemapped) == 0:
                print(f"No finemapped SNPs for locus {locname}, skipping")
                continue
            finemapped.to_csv(f"locus_{locname}.susie.finemapped", sep="\t", index=False)
        except:
            print(f"Error in {locname}")
        cmd = f'rm loci/locus_{locname}*'
        process = subprocess.Popen([cmd], close_fds=True, shell=True)
        process.wait(2.5)
        
if __name__ == "__main__":
    main()

