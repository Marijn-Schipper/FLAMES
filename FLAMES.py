from annotate import main as annotate
from optimize_FLAMES import main as optimize
from FLAMES_scoring import main as FLAMES_scoring
import os
import sys
import argparse
import pandas as pd

def splash_screen():
    print('\n*********************************************************************')
    print('*Fine-mapped Locus Assesment Model of Effector geneS (FLAMES)')
    print('* Version 1.1.2')
    print('* (C) 2023 Marijn Schipper')
    print('*********************************************************************')
    print()


def functional_annotation(args):
    parser = argparse.ArgumentParser(description="Annotate finemapped loci")
    parser.add_argument("-c", "--credsets_file", help="File containing credible set")
    parser.add_argument("-o", "--outdir", help="Output directory, if not specifying output filenames in an indexfile", required=False, default = None)
    parser.add_argument(
        "-l",
        "--GenomicRiskLoci",
        help="Path to the file that describes the boundaries of the locus that belongs to each inputted credible set. To score all genes in locus regardless of distance, also change -f flag (NOT RECOMMENDED).",
        required=False,
    )
    parser.add_argument(
        "-g",
        "--genes",
        help="Path to the file that describes the genes that belong to each locus. Needs to contain 'ensg' column and 'GenomicLocus' column, where the latter should match the locus in indexfile",
        required=False,
        default="BP mapping",
    )
    parser.add_argument(
        "-a",
        "--annotation_dir",
        help="Directory containing the annotation files",
        required=True,
    )
    parser.add_argument("-b", "--build", help="Genome build", default="GRCh37")
    parser.add_argument(
        "-p", "--pops", help="PoPS preds outputfile of correspondig GWAS", required=True
    )
    parser.add_argument(
        "-m",
        "--magma_z",
        help="Path to the file that describes the MAGMA gene level z-scores output",
        required=True,
    )
    parser.add_argument(
        "-mt",
        "--magma_tissue",
        help="Path to the file that describes the MAGMA tissue expression output",
        required=True,
    )
    parser.add_argument(
        "-tp",
        "--true_positives",
        help="Path to the file that describes the true positive genes in each locus",
        default = False
    )
    parser.add_argument(
        "-id",
        "--indexfile",
        help="File containing the locus no and its corresponing credible set file path",
    )
    parser.add_argument(
        "-f",
        "--filter",
        help="Drop all variants with distance > filter. Default is 750kb. WARNING: FLAMES is trained on 750kb, so precision estimate might be inacurate if increased",
        default=750000,
        type=int,
    )
    parser.add_argument(
        "-cv", "--cmd_vep", help="path to the vep executable", default=False
    )
    parser.add_argument(
        "-vc", "--vep_cache", help="path to the vep cache", default=False
    )
    parser.add_argument(
        "-sc",
        "--SNP_col",
        help="column name of SNP column in credset file",
        default="cred1",
    )
    parser.add_argument(
        "-pc",
        "--prob_col",
        help="column name of PIP column in credset file",
        default="prob1",
    )
    parser.add_argument(
        "-t", "--tabix", help="path to tabix if using local CADD scores", default=False
    )
    parser.add_argument(
        "-cf",
        "--CADD_file",
        help="path to CADD scores file if using local CADD scores, must match build of inputted credible variants",
        default=False,
    )
    parser.add_argument('-c95', '--credset_95', help='Input "FALSE" to not subset to 0.95 credible set', required=False, default=True)
    args = parser.parse_args(args)
    # At least one of the arguments is required
    if not (args.indexfile or args.credsets_file):
        parser.error("At least one of --credsets_file or --indexfile is required.")
    elif args.indexfile and args.credsets_file:
        parser.error("Only one of --credsets_file or --indexfile is allowed.")
    elif args.indexfile:
        if os.path.isfile(args.indexfile) == False:
            parser.error(
                f"The indexfile {args.indexfile} does not exist or is not a file."
            )
        elif args.outdir == None and not 'Annotfiles' in open(args.indexfile).readlines()[0]:
            parser.error(
                "When using an indexfile, an output directory must be specified or the column annotfiles must exist."
            )
    annotate(
        args.credsets_file,
        args.indexfile,
        args.annotation_dir,
        args.build,
        args.pops,
        args.magma_z,
        args.magma_tissue,
        args.outdir,
        args.SNP_col,
        args.prob_col,
        args.GenomicRiskLoci,
        args.genes,
        args.true_positives,
        args.filter,
        args.cmd_vep,
        args.vep_cache,
        args.tabix,
        args.CADD_file,
        args.credset_95,
    )
    return

def optimize_FLAMES(args):
    parser = argparse.ArgumentParser(
        description="optimize XGB with PoPS from annotated finemapped loci"
    )
    parser.add_argument(
        "-i",
        "--input_files",
        help="File containing the paths of all the inputfiles, can be the same files used to train your XGB model",
        required=True,
    )
    parser.add_argument(
        "-o", "--outdir", help="Output directory for figures and stats", required=True
    )
    parser.add_argument(
        "-m",
        "--modelpath",
        help="Path to trained model",
        required=False,
        default=os.path.join(os.path.dirname(__file__), "model"),
    )
    parser.add_argument(
        "-d",
        "--distance",
        help="Maximum inclusion distance of genes, default is 750kb, can be set to include all when set to 0",
        required=False,
        default=750000,
    )
    args = parser.parse_args(args)
    optimize(args.modelpath, args.input_files, args.outdir, args.distance)
    return


def FLAMES(args):
    parser = argparse.ArgumentParser(description="score finemapped loci with FLAMES")
    parser.add_argument(
        "-i",
        "--input_files",
        help="File containing the paths of all the inputfiles, or directory containing annotated files with annotated_ in the filename",
        required=False,
    )
    parser.add_argument('-id', '--indexfile', help='Tab/space delim. file containing the annotated loci input files under the column name Annotfiles ', required=False)
    parser.add_argument("-o", "--outdir", help="Output directory", required=True)
    parser.add_argument(
        "-f",
        "--filename",
        help="Output filename",
        required=False,
        default="FLAMES_scores",
    )
    parser.add_argument(
        "-d",
        "--distance",
        help="Maximum inclusion distance of genes, all annotated genes included when set to 0",
        required=False,
        default=750000,
    )
    parser.add_argument(
        "-w",
        "--weight",
        help="XGB weight to use, default is 0.725 XGBoost, 0.275 PoPS",
        required=False,
        default=0.725,
    )
    parser.add_argument(
        "-m",
        "--modelpath",
        help="Path to trained model",
        required=False,
        default=os.path.join(os.path.dirname(__file__), "model"),
    )
    args = parser.parse_args(args)
    if args.indexfile != None:
        if os.path.isfile(args.indexfile):
            try:
                inputfiles = list(pd.read_csv(args.indexfile, sep = "\t")['Annotfiles'])
            except:
                raise Exception("Indexfile not found or not in the right format")
        else:
            inputfiles = args.input_files
    elif args.input_files != None:
        try:
            inputfiles = open(args.input_files).read().splitlines()
        except:
            raise Exception("Inputfiles not found or not in the right format")
    FLAMES_scoring(
        args.modelpath,
        inputfiles,
        args.weight,
        args.distance,
        args.outdir,
        args.filename,
    )
    return


def run_MAGMA_tisue_type(args):
    return


def main():
    splash_screen()
    if len(sys.argv) < 2:
        raise Exception("State FLAMES function like so: python FLAMES.py annotate/FLAMES [args...]")
    command = sys.argv[1]
    args = sys.argv[2:]
    if command == "annotate":
        functional_annotation(args)
    elif command == "FLAMES":
        FLAMES(args)
    else:
        raise Exception(
            "Command not recognized. Please use annotate or FLAMES as your first argument."
        )
    return sys.exit(0)


if __name__ == "__main__":
    main()
