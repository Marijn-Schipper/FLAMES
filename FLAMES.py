from annotate import main as annotate
from train_FLAMES import main as train
from optimize_FLAMES import main as optimize
from FLAMES_scoring import main as FLAMES_scoring
import os
import sys
import argparse
import pandas as pd

def splash_screen():
    print('*********************************************************************')
    print('*Fine-mapped Locus Assesment Model of Effector geneS (FLAMES)')
    print('* Version 1.0.0')
    print('* (C) 2023 Marijn Schipper')
    print('*********************************************************************')
    print()


def functional_annotation(args):
    parser = argparse.ArgumentParser(description="Annotate finemapped loci")
    parser.add_argument("-c", "--credsets_file", help="File containing credible set")
    parser.add_argument("-o", "--outdir", help="Output directory", required=True)
    parser.add_argument(
        "-l",
        "--GenomicRiskLoci",
        help="Path to the file that describes the boundaries of the locus that belongs to each inputted credible set",
        required=False,
    )
    parser.add_argument(
        "-g",
        "--genes",
        help="Path to the file that describes the genes that belong to each locus",
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
        "-proc", "--processes", help="Number of processes to use", default=False
    )
    parser.add_argument(
        "-tp",
        "--true_positives",
        help="Path to the file that describes the true positive genes in each locus",
    )
    parser.add_argument(
        "-id",
        "--indexfile",
        help="File containing the locus no and its corresponing credible set file path",
    )
    parser.add_argument(
        "-f",
        "--filter",
        help="Filter the credible sets to only include those with a posterior probability of 0.95 or greater",
        action="store_true",
    )
    parser.add_argument(
        "-cv", "--cmd_vep", help="path to the vep executable", default=False
    )
    parser.add_argument(
        "-vc", "--vep_cache", help="path to the vep cache", default=False
    )
    parser.add_argument(
        "-vd", "--vep_docker", help="path to the vep docker", default=False
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
        args.vep_docker,
        args.tabix,
        args.CADD_file,
    )
    return


def train_xgb(args):
    parser = argparse.ArgumentParser(
        description="Train XGB model from annotated finemapped loci"
    )
    parser.add_argument(
        "-i",
        "--input_files",
        help="File containing the paths of all the inputfiles",
        required=True,
    )
    parser.add_argument(
        "-o", "--outdir", help="Output directory for the trained model", required=True
    )
    parser.add_argument(
        "-d",
        "--distance",
        help="Maximum inclusion distance of genes, default is 750kb, can be set to include all when set to 0",
        required=False,
        default=750000,
    )
    parser.add_argument(
        "-f",
        "--features",
        help="features used for training, default is all. RAW will only use raw features, rel will only use locus scaled features",
        required=False,
        default="all",
    )
    parser.add_argument(
        "-s",
        "--seed",
        help="random seed to use, default is 42, random will randomize the seed",
        required=False,
        default=42,
    )
    args = parser.parse_args(args)
    train(args.outdir, args.input_files, args.features, args.distance, args.seed)
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
        required=True,
    )
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
    FLAMES_scoring(
        args.modelpath,
        args.input_files,
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
        raise Exception("State FLAMES function like so: python FLAMES.py annotate/train/FLAMES/optimize [args...]")
    command = sys.argv[1]
    args = sys.argv[2:]
    if command == "annotate":
        functional_annotation(args)
    elif command == "train":
        train_xgb(args)
    elif command == "optimize":
        optimize_FLAMES(args)
    elif command == "FLAMES":
        FLAMES(args)
    else:
        raise Exception(
            "Command not recognized. Please use annotate, train, optimize or FLAMES as your first argument."
        )
    return sys.exit(0)


if __name__ == "__main__":
    main()
