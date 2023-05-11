import annotate_finemapped_fuma_filt as annotate
import os
import sys
import argparse
import pandas as pd

def create_loci_index_dict(credsets, indexfile=False):
    index = {}
    if indexfile == False:
        for credset in os.listdir(os.path.normpath(credsets)):
            if 'locus_' in credset:
                locus = int(credset.split('_')[-1].strip('.creds'))
            index[locus] = credset
    else:
        df = pd.read_csv(indexfile, delim_whitespace=True, header=None, names=['locus', 'credset'])
        index = dict(zip(df['locus'], df['credset']))
    return index


def functional_annotation(args):
    parser = argparse.ArgumentParser(description='Annotate finemapped loci')
    parser.add_argument('-c','--credsets_dir', help='Directory containing FUMA output')
    parser.add_argument('-o','--outdir', help='Output directory', required=True)
    parser.add_argument('-l','--GenomicRiskLoci', help='Path to the file that describes the boundaries of the locus that belongs to each inputted credible set', required=True)
    parser.add_argument('-g','--genes', help='Path to the file that describes the genes that belong to each locus', required=True)
    parser.add_argument('-tp', 'true_positives', help='Path to the file that describes the true positive genes in each locus')
    parser.add_argument('-id', '--indexfile', help='File containing the locus no and its corresponing credible set file path')
    parser.add_argument('-f', '--filter', help='Filter the credible sets to only include those with a posterior probability of 0.95 or greater', action='store_true')
    
    args = parser.parse_args()
    # At least one of the arguments is required
    if not (args.indexfile or args.credsets_dir):
        parser.error('At least one of --credsets_dir or --indexfile is required.')
    annotate.main(args.credsets_dir, args.outdir, args.GenomicRiskLoci, args.genes, args.true_positives, args.indexfile, args.filter)
    return

def MAGMA(args):
    return

def FLAMES(args):
    return



def main():
    if len(sys.argv) < 2:
        print("Usage: python FLAMES.py command [args...]")
        sys.exit(1)
    command = sys.argv[1]
    args = sys.argv[2:]
    if command == 'annotate':
        functional_annotation(args)
    elif command == 'MAGMA':
        MAGMA(args)
    elif command == 'FLAMES':
        FLAMES(args)
    else:
        print('Command not recognized. Please use annotate, MAGMA, or FLAMES as an argument.')


if __name__ == '__main__':
    main()

