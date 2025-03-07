# import modules
import pandas as pd

import os
from pyliftover import LiftOver
import numpy as np
from sys import argv
import sys
import time
import Query_api
import subprocess
import tempfile

# disable pandas helper messages
pd.options.mode.chained_assignment = None  # default='warn'


### Define helper functions
# Function to validate build and return a standardized string
def run_function(func, message, *args):
    try:
        func(*args)
    except:
        sys.exit(1, message)
    return True


def validate_build(build):
    if build.upper() == "HG19" or build.upper() == "GRCH37":
        build = "GRCH37"
    elif build.upper() == "HG38" or build.upper() == "GRCH38":
        build = "GRCH38"
    else:
        sys.exit("Invalid build. Please use build hg19 or hg38")
    return build


# Wrapper function that print function runtime
def runtime(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        func(*args, **kwargs)
        end = time.time()
        print(f"Runtime: {round(end - start, 2)}")

    return wrapper


# Annotate a credible set based on variant position within annotation range
def create_positional_mapping_matrix(df1, df2, colnames=["start", "end"], chr=False):
    if not chr == False:
        if chr == "X":
            chr = 23
        df1 = df1[df1["chr"] == chr]
    if "pos_37" in df2.columns:
        poscol = "pos_37"
    else:
        poscol = "pos"
    # Reshape and perform matrix factorization to check whether a SNP is within a genomic range
    starts = df1[colnames[0]].values.reshape(-1, 1)
    ends = df1[colnames[1]].values.reshape(-1, 1)
    positions = df2[poscol].values.reshape(1, -1)
    matches = np.logical_and(starts <= positions, ends >= positions)
    return matches


def liftover_row(chrom, pos, lo):
    newpos = lo.convert_coordinate(f"chr{chrom}", int(pos))
    try:
        newpos = newpos[0][1]
    except:
        newpos = np.nan
        print(f"\nUnable to liftover position: {chrom}:{pos}\n")
    return newpos


# Fuction to liftover a pandas column of genomic locations to different build
def liftover_df(lo, df, incols, outcol):
    df[outcol[0]] = df.apply(
        lambda row: liftover_row(f"{row[incols[0]]}", int(row[incols[1]]), lo),
        axis=1,
    )
    df = df.dropna(subset=[outcol[0]])
    if len(df) == 0:
        print("\nNo variants in credible set could be lifted over\n")
        return None
    return df


### Functions for annotation of credible sets
# Get the genes in a locus based on the FUMA genes file
def get_genes_in_locus(locus_no, path_to_genes_file, baseline_genes):
    genes = pd.read_csv(baseline_genes, sep="\t")
    if not "GenomicLocus" in genes.columns:
        gene_info = pd.read_csv(genes, sep="\t")
        relevant_genes = baseline_genes[
            baseline_genes["external_gene_name"].isin(genes["gene"])
        ]
        relevant_genes = pd.concat(
            [
                relevant_genes,
                baseline_genes[baseline_genes["symbol"].isin(genes["gene"])],
            ]
        )
        return relevant_genes
    genes["GenomicLocus"] = genes["GenomicLocus"].astype(str)
    genes["GenomicLocus"] = genes.GenomicLocus.str.split(":")
    genes = genes.explode("GenomicLocus")
    relevant_genes = pd.DataFrame()
    locusname = str(locus_no)
    genes_in_loc = genes[genes["GenomicLocus"] == locusname]
    relevant_genes = pd.concat([relevant_genes, genes_in_loc])
    if not "start" in relevant_genes.columns:
        gene_info = pd.read_csv(baseline_genes, sep="\t")
        relevant_genes = relevant_genes.merge(
            gene_info, left_on="gene", right_on="external_gene_name", how="left"
        )
        relevant_genes = relevant_genes.merge(
            gene_info, left_on="gene", right_on="symbol", how="left"
        )
    return relevant_genes


# Get the genes in locus based on physical location of the locus
def Genes_in_locus_bp_based(locus_no, GenomicRiskLoci, baseline_genes, lo=None):
    Locus = GenomicRiskLoci[
        GenomicRiskLoci["GenomicLocus"].astype(str) == str(locus_no)
    ]
    start = Locus["start"].iloc[0]
    end = Locus["end"].iloc[0]
    chrom = Locus["chr"].iloc[0]

    if chrom == "X":
        chrom = "23"
    # liftover to GRCh37 if build is GRCh38, retry slightly larger windows if a coordinate fails to convert
    if not lo == False:
        converted = False
        i = 0
        while not converted:
            new_start = lo.convert_coordinate(f"chr{chrom}", int(start) + i)
            if len(new_start) > 0:
                start = new_start[0][1]
                converted = True
            else:
                i += 50000
        converted = False
        i = 0
        while not converted:
            new_end = lo.convert_coordinate(f"chr{chrom}", int(end) - i)
            if len(new_end) > 0:
                end = new_end[0][1]
                converted = True
            else:
                i += 50000

    genes = pd.read_csv(baseline_genes, sep="\t", engine="pyarrow")
    genes["chromosome_name"] = genes["chromosome_name"].replace(["X", "x"], "23")

    genes = genes[
        (genes["end_position"] >= start)
        & (genes["start_position"] <= end)
        & (genes["chromosome_name"].astype(int) == int(chrom))
    ]
    genes = genes.rename(
        columns={
            "chromosome_name": "chr",
            "start_position": "start",
            "end_position": "end",
            "ensembl_gene_id": "ensg",
            "external_gene_name": "symbol",
        }
    )
    return genes


# Check if the genes are true positive training genes
def check_TP_genes(genes_in_locus, TP_genes):
    with open(TP_genes) as file:
        TP_genes = [line.rstrip() for line in file]
    ## check if genes in locus symbol or ensgid are in TP_genes and add column to df
    genes_in_locus = genes_in_locus[
        genes_in_locus["ensg"].isin(TP_genes) | genes_in_locus["symbol"].isin(TP_genes)
    ]
    return genes_in_locus


# Calculate the distances to the TSS and gene body from the credible set
def get_TSS_distances(creds, genes_in_locus, prob_col, build):
    distances = []
    weighted_distances = []
    TSS_distances = []
    weighted_TSS_distances = []
    poscol = "pos_37" if "pos_37" in creds.columns else "pos"
    centroid = sum(creds[poscol] * creds[prob_col]) / sum(creds[prob_col])
    most_likely_snp = creds.loc[creds[prob_col].idxmax()][poscol]
    for index, row in genes_in_locus.iterrows():
        if row["strand"] == "+":
            TSS_distances.append(abs(row["start"] - most_likely_snp))
            weighted_TSS_distances.append(abs(row["start"] - centroid))
        else:
            TSS_distances.append(abs(most_likely_snp - row["end"]))
            weighted_TSS_distances.append(abs(centroid - row["end"]))
        if row["start"] < centroid < row["end"]:
            distances.append(0)
            weighted_distances.append(0)
        else:
            distances.append(
                min(
                    abs(row["start"] - most_likely_snp),
                    abs(row["end"] - most_likely_snp),
                )
            )
            weighted_distances.append(
                min(abs(row["start"] - centroid), abs(row["end"] - centroid))
            )
    genes_in_locus["TSS_distance"] = TSS_distances
    genes_in_locus["weighted_TSS_distance"] = weighted_TSS_distances
    genes_in_locus["distance"] = distances
    genes_in_locus["weighted_distance"] = weighted_distances
    return genes_in_locus


# Create a dictionary of mamga tissue type relevance for GTEx eQTLs
def create_relevance_dict(magma_scores):
    magma_scores["score"] = magma_scores["P"].rank(ascending=False) / len(
        magma_scores["P"]
    )
    threshold = 0.05 / len(magma_scores["score"])
    magma_scores["score"] = magma_scores.apply(
        lambda row: row["score"] * 0.5
        if row["P"] > 0.05
        else row["score"] * 2
        if row["P"] <= threshold
        else row["score"],
        axis=1,
    )
    magma_scores["VARIABLE"] = magma_scores["VARIABLE"].str.slice(stop=27)
    magma_scores_dict = dict(zip(magma_scores["VARIABLE"], magma_scores["score"]))
    return magma_scores_dict


# Get the GTEx eQTLs in a locus
def get_GTEX_eQTLs(creds, genes, prob_col, finemapped_gtex, build, magma_scores=False):
    # prep data for input build
    if build.upper() == "GRCH38":
        creds["eqtl_id"] = creds.apply(
            lambda row: f"chr{row['chr']}_{row['pos']}_{row['a1']}_{row['a2']}_b38",
            axis=1,
        )
        creds["eqtl_id2"] = creds.apply(
            lambda row: f"chr{row['chr']}_{row['pos']}_{row['a2']}_{row['a1']}_b38",
            axis=1,
        )
        header = "eQTL"
    else:
        creds["eqtl_id"] = creds.apply(
            lambda row: f"{row['chr']}_{row['pos']}_{row['a1']}_{row['a2']}_b37", axis=1
        )
        creds["eqtl_id2"] = creds.apply(
            lambda row: f"{row['chr']}_{row['pos']}_{row['a2']}_{row['a1']}_b37", axis=1
        )
        header = "eQTL_hg37"
    # GTEX get CLPP
    CLPP_sum = []
    CLPP_max = []
    tissue_specific_CLPP_sum = []
    tissue_specific_CLPP_max = []
    for index, row in genes.iterrows():
        maxval, sumval = 0, 0
        genepath = os.path.join(finemapped_gtex, row["ensg"].split(".")[0] + ".parquet")
        if not os.path.isfile(genepath):
            CLPP_sum.append(sumval)
            CLPP_max.append(maxval)
            tissue_specific_CLPP_sum.append(sumval)
            tissue_specific_CLPP_max.append(maxval)
            continue
        eqtls = pd.read_parquet(os.path.join(genepath), engine="fastparquet")

        overlap = pd.merge(
            creds, eqtls, left_on="eqtl_id", right_on=header, how="inner"
        )
        overlap = pd.concat(
            [
                overlap,
                pd.merge(
                    creds, eqtls, left_on="eqtl_id2", right_on=header, how="inner"
                ),
            ]
        )
        if overlap.shape[0] > 0:
            sumval = sum(list(overlap["Probability"] * overlap[prob_col]))
            maxval = max(list(overlap["Probability"] * overlap[prob_col]))
        CLPP_sum.append(sumval)
        CLPP_max.append(maxval)
        if type(magma_scores) == dict and overlap.shape[0] > 0:
            overlap["TISSUE"] = overlap["TISSUE"].str.slice(stop=27)
            overlap["scores"] = overlap["TISSUE"].map(magma_scores)
            overlap["Probability"] = overlap["Probability"] * overlap["scores"]
            sumval = sum(list(overlap["Probability"] * overlap[prob_col]))
            maxval = max(list(overlap["Probability"] * overlap[prob_col]))
            tissue_specific_CLPP_sum.append(sumval)
            tissue_specific_CLPP_max.append(maxval)
        else:
            tissue_specific_CLPP_sum.append(sumval)
            tissue_specific_CLPP_max.append(maxval)
    genes["CLPP_GTEX_eQTL_sum"] = CLPP_sum
    genes["CLPP_GTEX_eQTL_max"] = CLPP_max
    if type(magma_scores) == dict:
        genes["CLPP_GTEX_tissue_weighted_eQTL_sum"] = tissue_specific_CLPP_sum
        genes["CLPP_GTEX_tissue_weighted_eQTL_max"] = tissue_specific_CLPP_max
    return genes


# Get the GTEx eQTLs in a locus
def get_eqtlgen_eQTLs(creds, genes, prob_col, eQTL_dir, build):
    chrom = creds.iloc[0]["chr"]
    if build.upper() == "GRCH38":
        creds["eqtl_id"] = creds.apply(
            lambda row: f"chr{row['chr']}_{row['pos']}_{row['a1']}_{row['a2']}_b38",
            axis=1,
        )
        creds["eqtl_id2"] = creds.apply(
            lambda row: f"chr{row['chr']}_{row['pos']}_{row['a2']}_{row['a1']}_b38",
            axis=1,
        )
        creds["eqtl_id3"] = creds.apply(
            lambda row: f"{row['chr']}_{row['pos']}_X_X_b38", axis=1
        )
        header = "hg38"
    else:
        creds["eqtl_id"] = creds.apply(
            lambda row: f"{row['chr']}_{row['pos']}_{row['a1']}_{row['a2']}_b37", axis=1
        )
        creds["eqtl_id2"] = creds.apply(
            lambda row: f"{row['chr']}_{row['pos']}_{row['a2']}_{row['a1']}_b37", axis=1
        )
        creds["eqtl_id3"] = creds.apply(
            lambda row: f"{row['chr']}_{row['pos']}_X_X_b37", axis=1
        )
        header = "hg37"
    CLPP_max = []
    CLPP_sum = []
    for index, row in genes.iterrows():
        maxval, sumval = 0, 0
        genepath = os.path.join(
            f"{eQTL_dir}", f"chr{chrom}", f"{row['ensg'].upper()}.parquet"
        )
        if not os.path.isfile(genepath):
            CLPP_sum.append(sumval)
            CLPP_max.append(maxval)
            continue
        eqtls = pd.read_parquet(os.path.join(genepath), engine="fastparquet")
        overlap = pd.merge(
            creds, eqtls, left_on="eqtl_id", right_on=header, how="inner"
        )
        overlap = pd.concat(
            [
                overlap,
                pd.merge(
                    creds, eqtls, left_on="eqtl_id2", right_on=header, how="inner"
                ),
            ]
        )
        overlap = pd.concat(
            [
                overlap,
                pd.merge(
                    creds, eqtls, left_on="eqtl_id3", right_on=header, how="inner"
                ),
            ]
        )
        if overlap.shape[0] > 0:
            sumval = sum(list(overlap["Causal_Post._Prob."] * overlap[prob_col]))
            maxval = max(list(overlap["Causal_Post._Prob."] * overlap[prob_col]))
        CLPP_sum.append(sumval)
        CLPP_max.append(maxval)
    genes["CLPP_eQTLgen_sum"] = CLPP_sum
    genes["CLPP_eQTLgen_max"] = CLPP_max
    return genes


# get genes from eQTL catalog
def get_QTL_catalog(creds, genes, prob_col, eQTL_catalog, build):
    if build.upper() == "GRCH38":
        creds["eqtl_id"] = creds.apply(
            lambda row: f"chr{row['chr']}_{row['pos']}_{row['a1']}_{row['a2']}",
            axis=1,
        )
        creds["eqtl_id2"] = creds.apply(
            lambda row: f"chr{row['chr']}_{row['pos']}_{row['a2']}_{row['a1']}",
            axis=1,
        )
        header = "QTL_b38"
    else:
        creds["eqtl_id"] = creds.apply(
            lambda row: f"chr{row['chr']}_{row['pos']}_{row['a1']}_{row['a2']}", axis=1
        )
        creds["eqtl_id2"] = creds.apply(
            lambda row: f"chr{row['chr']}_{row['pos']}_{row['a2']}_{row['a1']}", axis=1
        )
        header = "QTL_b37"
    for qtl in ["eQTL", "rQTL"]:
        QTL_max_CLPP_max = []
        QTL_max_CLPP_sum = []
        QTL_mean_CLPP_max = []
        QTL_mean_CLPP_sum = []
        for index, row in genes.iterrows():
            maxval, sumval, max2val, sum2val = 0, 0, 0, 0
            genepath = os.path.join(
                f"{os.path.join(eQTL_catalog, qtl)}", f"{row['ensg'].upper()}.parquet"
            )
            if os.path.isfile(genepath):
                eqtls = pd.read_parquet(genepath, engine="fastparquet")
                overlap = pd.merge(
                    creds, eqtls, left_on="eqtl_id", right_on=header, how="inner"
                )
                overlap = pd.concat(
                    [
                        overlap,
                        pd.merge(
                            creds,
                            eqtls,
                            left_on="eqtl_id2",
                            right_on=header,
                            how="inner",
                        ),
                    ]
                )
                if overlap.shape[0] > 0:
                    sumval = sum(list(overlap["max"] * overlap[prob_col]))
                    maxval = max(list(overlap["max"] * overlap[prob_col]))
                    max2val = max(list(overlap["mean"] * overlap[prob_col]))
                    sum2val = sum(list(overlap["mean"] * overlap[prob_col]))
            QTL_max_CLPP_max.append(maxval)
            QTL_max_CLPP_sum.append(sumval)
            QTL_mean_CLPP_max.append(max2val)
            QTL_mean_CLPP_sum.append(sum2val)
        genes[f"max_CLPP_{qtl}_eQTLCatalog_sum"] = QTL_max_CLPP_sum
        genes[f"max_CLPP_{qtl}_eQTLCatalog_max"] = QTL_max_CLPP_max
        genes[f"mean_CLPP_{qtl}_eQTLCatalog_sum"] = QTL_mean_CLPP_sum
        genes[f"mean_CLPP_{qtl}_eQTLCatalog_max"] = QTL_mean_CLPP_max
    return genes


def get_Promoters(prom_dir, genes, creds, prob_col, build):
    sums = []
    maxs = []
    proms = pd.read_parquet(
        os.path.join(prom_dir, f"Proms_per_transcript_{build}.parquet")
    )
    for gene in genes["ensg"]:
        prom_df = proms[proms["ensg"] == gene]
        if len(prom_df) == 0:
            sums.append(0)
            maxs.append(0)
            continue
        ## create a matrix that maps which variant is in which captured region
        matrix = create_positional_mapping_matrix(prom_df, creds)
        matrix2 = pd.DataFrame(np.outer(np.ones(len(prom_df)), creds[prob_col]))
        matrix2 = matrix2.where(matrix, other=0)
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["Promoter_sum"] = sums
    genes["Promoter_max"] = maxs
    return genes


# Get the GTEx eQTLs in a locus
def get_VEP(creds, genes, prob_col, build):
    VEP_dict = {"HIGH": 1, "MODERATE": 0.6, "LOW": 0.4, "MODIFIER": 0.1}
    VEPs = creds.apply(
        lambda row: Query_api.query_VEP(
            row["chr"], row["pos"], row["a1"], row["a2"], build
        ),
        axis=1,
    )
    if any(VEPs.apply(lambda x: isinstance(x, str))):
        return "Cancel annotation due to timeout in VEP"
    vars = []
    found_genes = []
    consequences = []
    for i, query in enumerate(VEPs):
        variants = query[0]
        PiP = creds[prob_col][i]
        if not "transcript_consequences" in variants.keys():
            continue
        for variant in variants["transcript_consequences"]:
            if not "gene_id" in variant.keys():
                continue
            consequence = float(PiP) * float(VEP_dict[variant["impact"]])
            vars.append(i)
            found_genes.append(variant["gene_id"])
            consequences.append(consequence)
    vep_df = pd.DataFrame()
    vep_df["ensg"] = found_genes
    vep_df["consequence"] = consequences
    vep_df['variant'] = vars  
    vep_df = vep_df.groupby(["ensg", "variant"])["consequence"].max().reset_index()  
    vep_df["VEP_max"] = vep_df.groupby("ensg")["consequence"].transform("max")
    vep_df["VEP_sum"] = vep_df.groupby("ensg")["consequence"].transform("sum")
    vep_df = vep_df[["ensg", "VEP_sum", "VEP_max"]]
    genes = genes.merge(vep_df, on="ensg", how="left")
    genes.fillna(0, inplace=True)
    return genes


# run vep in environment
def run_vep_within_environment(
    input_file, output_file, vep_path, vep_cache, build,
):
    if build.upper() == "GRCH38" or build.upper() == "HG38":
        build = "GRCh38"
    else:
        build = "GRCh37"
    vep_command = f"{vep_path} -i {input_file} -o {output_file} --cache {vep_cache} --assembly {build} --offline --force_overwrite"
    subprocess.run(vep_command, shell=True, check=True)
    return


def cmd_VEP(
    creds, genes, prob_col, build, VEP_path, VEP_cache, outdir
):
    tmpname = tempfile.mktemp()
    tmp_file_path = os.path.join(f"{tmpname}.vcf")
    output_file_path = os.path.join(f"{tmpname}_output.vcf")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    output_lines = creds.copy()
    output_lines = output_lines.sort_values(["chr", "pos"])
    output_lines[['a1', 'a2']] = output_lines[['a1', 'a2']]


    output_lines = output_lines.apply(
        lambda row: [
            f"{row['chr']} {row['pos']} {row['pos']} {row['a1']}/{row['a2']} +",
            f"{row['chr']} {row['pos']} {row['pos']} {row['a1']}/{row['a2']} -",
        ],
        axis=1,
    )
    # Flatten the list of output lines
    output_lines = [line for variant_lines in output_lines for line in variant_lines]

    # Write the lines to a file
    with open(tmp_file_path, "w") as file:
        file.write("\n".join(output_lines))

    run_vep_within_environment(
        tmp_file_path, output_file_path, VEP_path, VEP_cache, build,
    )
    data = []
    with open(output_file_path, "r") as file:
        for line in file:
            if not line.startswith("##"):
                data.append(line.strip().split("\t"))
    
    vep_df = pd.DataFrame(data[1:], columns=data[0])
    vep_df["Extra"] = vep_df["Extra"].str.extract(r"IMPACT=(\w+);")
    vep_df["Extra"] = vep_df["Extra"].fillna("None")
    VEP_dict = {"HIGH": 1, "MODERATE": 0.6, "LOW": 0.4, "MODIFIER": 0.1, "None": 0}
    vep_df["VEP"] = vep_df["Extra"].apply(lambda x: VEP_dict[x])
    vep_df = vep_df.groupby(["Gene", "#Uploaded_variation"])["VEP"].max().reset_index()
    creds["VEP_merge_id"] = creds.apply(
        lambda row: f"{row['chr']}_{row['pos']}_{row['a1']}/{row['a2']}",
        axis=1,
    )
    vep_df = vep_df.merge(
        creds, left_on="#Uploaded_variation", right_on="VEP_merge_id", how="left"
    )
    creds.drop(columns=["VEP_merge_id"], inplace=True)
    vep_df["VEP_weighted"] = vep_df["VEP"] * vep_df[prob_col]
    vep_df["VEP_max"] = vep_df.groupby("Gene")["VEP_weighted"].transform("max")
    vep_df["VEP_sum"] = vep_df.groupby("Gene")["VEP_weighted"].transform("sum")
    vep_df = vep_df[["Gene", "VEP_sum", "VEP_max"]]
    genes = genes.merge(vep_df, left_on="ensg", right_on="Gene", how="left")
    genes.drop(columns=["Gene"], inplace=True)
    genes.fillna(0, inplace=True)
    for f in [
        tmp_file_path,
        output_file_path,
        output_file_path + "_summary.html",
        output_file_path + "_warnings.txt",
    ]:
        if os.path.exists(f):
            os.remove(f)
    return genes


# Get the CADD scores within each gene in the locus
def get_CADD(creds, genes, prob_col, build):
    creds["CADD"] = creds.apply(
        lambda row: Query_api.query_CADD(
            row["chr"], row["pos"], row["a1"], row["a2"], build
        ),
        axis=1,
    )
    matrix = create_positional_mapping_matrix(genes, creds)
    creds["CADD_weighted"] = creds["CADD"] * creds[prob_col]
    genes["CADD_sum"] = np.dot(matrix, creds["CADD_weighted"])
    max_values = np.nanmax(np.where(matrix, creds["CADD_weighted"], 0), axis=1)
    genes["CADD_max"] = max_values
    return genes


def tabix_CADD(chr, pos, a1, a2, CADD_path, tabix_path):
    region_of_interest = f"{chr}:{pos}-{pos}"

    # Run the tabix command and capture the output
    command = f"{tabix_path} {CADD_path} {region_of_interest}"
    result = subprocess.run(command, capture_output=True, shell=True)
    tabix_output = result.stdout.decode("utf-8")
    header = ["chr", "pos", "a1", "a2", "raw", "PHRED"]
    tabix_output = tabix_output.split("\n")
    tabix_output = [x.split("\t") for x in tabix_output if len(x.strip()) > 0]
    df = pd.DataFrame(tabix_output, columns=header)  # ,dtype={'p':'string'})
    if len(df) == 0:
        return 0
    ndf = df[(df["a1"] == a1) & (df["a2"] == a2)]
    if len(ndf) == 1:
        return float(ndf["PHRED"].iloc[0])
    ndf = df[(df["a1"] == a2) & (df["a2"] == a1)]
    if len(ndf) == 1:
        return float(ndf["PHRED"].iloc[0])
    return 0


def get_CADD_cmd(creds, genes, prob_col, CADD_path, tabix_path):
    creds["CADD"] = creds.apply(
        lambda row: tabix_CADD(
            row["chr"], row["pos"], row["a1"], row["a2"], CADD_path, tabix_path
        ),
        axis=1,
    )
    matrix = create_positional_mapping_matrix(genes, creds)
    creds["CADD_weighted"] = creds["CADD"] * creds[prob_col]
    genes["CADD_sum"] = np.dot(matrix, creds["CADD_weighted"])
    max_values = np.nanmax(np.where(matrix, creds["CADD_weighted"], 0), axis=1)
    genes["CADD_max"] = max_values
    return


# Get the Jung HiC interaction scores for each gene in the locus
def get_jung_HiC(jung, genes, creds, prob_col):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(jung, gene + ".parquet")):
            sums.append(0)
            maxs.append(0)
            continue
        jungs = pd.read_parquet(
            os.path.join(jung, gene + ".parquet"), engine="fastparquet"
        )
        ## split jungs 'raw' on '.' and split into three columns called chr, start and end
        jungs[["tmp", "start", "end"]] = jungs["target"].str.split(".", expand=True)
        jungs["p"] = jungs["dist_pvalue"]
        jungs["start"] = jungs["start"].astype(int)
        jungs["end"] = jungs["end"].astype(int)
        jungs = jungs.drop(jungs[jungs["p"] >= 1000].index)
        ## create a matrix that maps which variant is in which captured region
        matrix = create_positional_mapping_matrix(jungs, creds)
        matrix2 = pd.DataFrame(np.outer(jungs["p"], creds[prob_col]))
        matrix2 = matrix2.where(matrix, other=0)
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["Jung_HiC_sum"] = sums
    genes["Jung_HiC_max"] = maxs
    return genes


# Get the Javierre HiC interaction scores for each gene in the locus
def get_javierre_HiC(javierre, genes, creds, prob_col):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(javierre, gene + ".parquet")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        javierres = pd.read_parquet(
            os.path.join(javierre, gene + ".parquet"), engine="fastparquet"
        )
        celltypes = [
            "Mon",
            "Mac0",
            "Mac1",
            "Mac2",
            "Neu",
            "MK",
            "EP",
            "Ery",
            "FoeT",
            "nCD4",
            "tCD4",
            "aCD4",
            "naCD4",
            "nCD8",
            "tCD8",
            "nB",
            "tB",
        ]
        ## create a max and sum column in javierres based on the celltypes values
        javierres["max"] = javierres[celltypes].max(axis=1)
        javierres["sum"] = javierres[celltypes].sum(axis=1)
        # melt the dataframe to have each celltype measurement into individual rows
        matrix = create_positional_mapping_matrix(
            javierres, creds, ["oeStart", "oeEnd"]
        )
        matrix2 = pd.DataFrame(np.outer(javierres["sum"], creds[prob_col])).where(
            matrix, other=0
        )
        matrix3 = pd.DataFrame(np.outer(javierres["max"], creds[prob_col])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix3.values.max())
    genes["Javierre_HiC_sum"] = sums
    genes["Javierre_HiC_max"] = maxs
    return genes


# Get the ABC EP interaction scores for each gene in the locus
def get_HACER(HACER, genes, creds, prob_col):
    for method in ["CAGE", "PRO-seq_GRO-seq"]:
        sums = []
        maxs = []
        for gene in genes["ensg"]:
            HACERs = pd.DataFrame()
            for d in ["F5", "closest", "4D"]:
                try:
                    tmpdf = pd.read_parquet(
                        os.path.join(HACER, method, d, gene + ".parquet"),
                        engine="fastparquet",
                    )
                    HACERs = pd.concat([HACERs, tmpdf])
                except:
                    pass
            if len(HACERs) == 0:
                sums.append(0.0)
                maxs.append(0.0)
                continue
            HACERs["chr"] = HACERs["chr"].str.replace("chr", "")
            if not "density" in HACERs.columns:
                HACERs["density"] = 1
            HACERs["start"] = HACERs["start"].astype(int)
            HACERs["end"] = HACERs["end"].astype(int)
            matrix = create_positional_mapping_matrix(HACERs, creds)
            matrix2 = pd.DataFrame(np.outer(HACERs["density"], creds[prob_col])).where(
                matrix, other=0
            )
            sums.append(matrix2.values.sum())
            maxs.append(matrix2.values.max())
        genes[f"HACER_{method}_sum"] = sums
        genes[f"HACER_{method}_max"] = maxs

    return genes


# Get the ABC EP interaction scores for each gene in the locus
def get_CICERO(cicero, genes, creds, prob_col):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(cicero, gene + ".parquet")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        cics = pd.read_parquet(
            os.path.join(cicero, gene + ".parquet"), engine="fastparquet"
        )
        cics["chr"] = cics["chr"].str.replace("chr", "")
        matrix = create_positional_mapping_matrix(cics, creds)
        matrix2 = pd.DataFrame(np.outer(cics["cor"], creds[prob_col])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes[f"cicero_sum"] = sums
    genes[f"cicero_max"] = maxs
    return genes


# Get the ABC EP interaction scores for each gene in the locus
def get_ABC_EP(ABC_EP, genes, creds, prob_col):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(ABC_EP, gene + ".parquet")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        ABCs = pd.read_parquet(
            os.path.join(ABC_EP, gene + ".parquet"), engine="fastparquet"
        )
        matrix = create_positional_mapping_matrix(ABCs, creds)
        matrix2 = pd.DataFrame(np.outer(ABCs["ABC.Score"], creds[prob_col])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["ABC_EP_sum"] = sums
    genes["ABC_EP_max"] = maxs
    return genes


# Get the ABC CRISPR interaction scores for each gene in the locus
def get_ABC_CRISPR(ABC_CRISPR, genes, creds, prob_col):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(ABC_CRISPR, gene + ".parquet")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        ABCs = pd.read_parquet(
            os.path.join(ABC_CRISPR, gene + ".parquet"), engine="fastparquet"
        )
        matrix = create_positional_mapping_matrix(ABCs, creds)
        matrix2 = pd.DataFrame(np.outer(ABCs["ABC_max"], creds[prob_col])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["ABC_CRISPR_sum"] = sums
    genes["ABC_CRISPR_max"] = maxs
    return genes


# Get the FANTHOM5 interaction scores for each gene in the locus
def get_F5(F5_ep, genes, creds, prob_col):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(F5_ep, gene + ".parquet")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        F5s = pd.read_parquet(
            os.path.join(F5_ep, gene + ".parquet"), engine="fastparquet"
        )
        F5s[["chr", "start", "end"]] = F5s["enhancer"].str.split(r"[:|-]", expand=True)
        F5s["start"] = F5s["start"].astype(int)
        F5s["end"] = F5s["end"].astype(int)
        matrix = create_positional_mapping_matrix(F5s, creds)
        matrix2 = pd.DataFrame(np.outer(F5s["correlation"], creds[prob_col])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["F5_ep_sum"] = sums
    genes["F5_ep_max"] = maxs
    return genes


# Get Genehancer predictions for each gene in locus
def get_Genehancer(genehancer, genes, creds, prob_col):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(genehancer, gene + ".parquet")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        genehancer_scores = pd.read_parquet(
            os.path.join(genehancer, gene + ".parquet"), engine="fastparquet"
        )
        matrix = create_positional_mapping_matrix(genehancer_scores, creds)
        matrix2 = pd.DataFrame(
            np.outer(genehancer_scores["score"], creds[prob_col])
        ).where(matrix, other=0)
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["GeneHancer_max"] = maxs
    genes["GeneHancer_sum"] = sums
    return genes


# Get Roadmap Epigenomics scores for each gene in locus
def get_RoadmapEpi(Roadmap, genes, creds, prob_col):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(Roadmap, gene + ".parquet")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        Epi = pd.read_parquet(
            os.path.join(Roadmap, gene + ".parquet"), engine="fastparquet"
        )
        matrix = create_positional_mapping_matrix(Epi, creds)
        matrix2 = pd.DataFrame(np.outer(Epi["max"], creds[prob_col])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["RoadmapEpi_sum"] = sums
    genes["RoadmapEpi_max"] = maxs
    return genes


# Get EpiMap scores for each gene in locus
def get_EpiMap(EpiMap, genes, creds, prob_col):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(EpiMap, gene + ".parquet")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        Epi = pd.read_parquet(
            os.path.join(EpiMap, gene + ".parquet"), engine="fastparquet"
        )
        matrix = create_positional_mapping_matrix(Epi, creds)
        matrix2 = pd.DataFrame(np.outer(Epi["corr"], creds[prob_col])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["EpiMap_sum"] = sums
    genes["EpiMap_max"] = maxs
    return genes


# Run all annotation funcs for each annotation dataset
def annotate_credset(
    genes,
    creds,
    build,
    Ann_path,
    cmd_vep,
    vep_cache,
    tabix,
    CADD_file,
    outdir,
    prob_col,
    fname,
    magma=False,
    trainset=False,
    filter=750000,
):
    genes = get_TSS_distances(creds, genes, prob_col, build)
    if trainset:
        if not check_dist_to_causal(genes, dist=750000):
            exit("The true positive causal gene is too far from lead SNP for {fname}")
    genes = genes[genes["weighted_distance"] <= filter]
    if len(genes) == 0:
        print(f"\nNo genes after filtering SNPs on distance, exiting\n")
        return genes
    start_time = time.time()
    if cmd_vep == False:
        genes = get_VEP(creds, genes, prob_col, build)
        if len(genes) == 0:
            print("\nError in querying VEP, exiting\n")
            return genes
    else:
        genes = cmd_VEP(
            creds, genes, prob_col, build, cmd_vep, vep_cache, outdir
        )
        if len(genes) == 0:
            print("\nError in querying VEP, this could be due to the VEP API server. Please try again later or contact our GitHub, exiting\n")
            return genes
    if type(genes) == str:
        print("\nError in querying VEP, exiting\n")
        return genes
    if CADD_file == False or tabix == False:
        genes = get_CADD(creds, genes, prob_col, build)
        if len(genes) == 0:
            print("\nError in querying CADD, this could be due to the VEP API server. Please try again later or contact our GitHub, exiting\n")
            return genes
    else:
        get_CADD_cmd(creds, genes, prob_col, CADD_file, tabix)
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_eqtlgen_eQTLs(
        creds, genes, prob_col, os.path.join(Ann_path, "eQTLGen"), build
    )
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_GTEX_eQTLs(
        creds,
        genes,
        prob_col,
        os.path.join(Ann_path, "GTEx_v8_finemapping_CAVIAR"),
        build,
        magma,
    )
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_QTL_catalog(
        creds, genes, prob_col, os.path.join(Ann_path, "eQTL_catalog"), build
    )
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_HACER(os.path.join(Ann_path, "HACER"), genes, creds, prob_col)
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_CICERO(
        os.path.join(Ann_path, "Cicero_whole_blood"), genes, creds, prob_col
    )
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_jung_HiC(os.path.join(Ann_path, "Jung_PCHiC"), genes, creds, prob_col)
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_javierre_HiC(
        os.path.join(Ann_path, "Javierre_PCHiC"), genes, creds, prob_col
    )
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_ABC_EP(os.path.join(Ann_path, "ABC_EP"), genes, creds, prob_col)
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_ABC_CRISPR(os.path.join(Ann_path, "ABC_CRISPR"), genes, creds, prob_col)
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_F5(os.path.join(Ann_path, "EP_correlation_F5"), genes, creds, prob_col)
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_Genehancer(os.path.join(Ann_path, "GeneHancer"), genes, creds, prob_col)
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_RoadmapEpi(os.path.join(Ann_path, "RoadMapEpi"), genes, creds, prob_col)
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_EpiMap(os.path.join(Ann_path, "EpiMap"), genes, creds, prob_col)
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes = get_Promoters(
        os.path.join(Ann_path, "Promoter_regions"), genes, creds, prob_col, build
    )
    return genes


# Get PoPS and MAGMA scores for each gene in locus
def POPS_MAGMA_annotation(genes, magma, PoPS):
    magma = magma[["GENE", "ZSTAT"]]
    genes = pd.merge(genes, magma, left_on="ensg", right_on="GENE", how="left")
    genes = genes.drop(columns=["GENE"])
    genes = genes.rename(columns={"ZSTAT": "MAGMA_Z"}).fillna(0)
    PoPS = PoPS[["ENSGID", "PoPS_Score"]]
    genes = pd.merge(genes, PoPS, left_on="ensg", right_on="ENSGID", how="left")
    genes = genes.drop(columns=["ENSGID"])
    return genes


# Transform all scores to locus relative scores
def create_relative_annotation(genes):
    for feature in genes.columns:
        if len(genes) < 1:
            continue
        if "rel_" in feature:
            continue
        if "distance" in feature:
            maxval = max(genes[feature].astype(float))
            minval = min(genes[feature].astype(float))
            if maxval == minval:
                genes["rel_" + feature] = 1
            else:
                genes["rel_" + feature] = maxval - genes[feature]
                genes["rel_" + feature] = genes["rel_" + feature] / float(
                    maxval - minval
                )
        elif any(x in feature for x in ["max", "sum", "PoPS", "MAGMA_Z"]):
            if feature == "sum":
                continue
            maxval = max(genes[feature].astype(float))
            minval = min(genes[feature].astype(float))
            if maxval == minval and maxval > 0.0:
                genes["rel_" + feature] = 1
            elif max(genes[feature]) == 0.0:
                genes[f"rel_{feature}"] = genes[feature].copy()
            else:
                genes[f"rel_{feature}"] = genes[feature].astype(float) - minval
                genes[f"rel_{feature}"] = genes[f"rel_{feature}"].astype(float) / float(
                    maxval - minval
                )
    return genes


# check if loci definitions are provided
def check_risk_loci(GenomicRiskLoci, creds, locno, filter):
    if not isinstance(GenomicRiskLoci, type(pd.DataFrame())):
        GenomicRiskLoci = pd.DataFrame()
        GenomicRiskLoci["GenomicLocus"] = [locno]
        GenomicRiskLoci["chr"] = min(creds["chr"])
        GenomicRiskLoci["start"] = max([0, min(creds["pos"]) - filter])
        GenomicRiskLoci["end"] = max(creds["pos"]) + filter
    return GenomicRiskLoci


# generate 95% credible set
def create_95perc_credset(creds, prob_col):
    sorted_data = creds.sort_values(prob_col)
    if sum(sorted_data[prob_col]) < 0.95:
        print("\nWARNING: The sum of the probabilities is less than 0.95")
        return(creds)
    if sum(sorted_data[prob_col]) > 1.0:
        print("\nWARNING: The sum of the probabilities is greater than 1.0, rescaling")
        sorted_data[prob_col] = sorted_data[prob_col] / sum(sorted_data[prob_col])
    sorted_data["cumulative_sum"] = sorted_data[prob_col].cumsum()
    creds = sorted_data[sorted_data["cumulative_sum"] >= 0.05]
    creds.reset_index(drop=True, inplace=True)
    creds = creds.drop(columns=["cumulative_sum"])
    return creds


# annotate locus metadata
def create_locus_metadata(genes, creds, prob_col):
    genes["genes_in_locus"] = len(genes)
    genes["credset_size"] = len(creds)
    genes["genes_within_50kb"] = len(genes[genes["distance"] <= 50000])
    genes["genes_within_100kb"] = len(genes[genes["distance"] <= 100000])
    genes["genes_within_250kb"] = len(genes[genes["distance"] <= 250000])
    genes["Highest_PIP"] = creds[prob_col].max()
    return genes


# Process entire locus
def full_annotation_of_credset(
    path_to_credset,
    build,
    GenomicRiskLoci,
    Ann_path,
    outdir,
    lo,
    magma_z,
    PoPS,
    genes,
    SNP_col="cred1",
    prob_col="prob1",
    magma_scores=False,
    trainset=False,
    filter=False,
    cmd_vep=False,
    vep_cache=False,
    tabix=False,
    CADD_file=False,
    c95=True
):
    # load credset
    locno = path_to_credset[1]
    given_outname = path_to_credset[2]
    path_to_credset = path_to_credset[0]
    if given_outname == None:
        outfile = os.path.join(
            outdir, f"FLAMES_annotated_{os.path.basename(path_to_credset)}"
        )
    else:
        outfile = given_outname
        try:
            outdir = os.path.dirname(outfile)
        except:
            outdir = os.getcwd()
    if outdir == "":
        outdir = os.getcwd()
    if os.path.exists(outfile):
        print(f"\nAnnotation file {outfile} already exists")
        return
    creds = pd.read_csv(path_to_credset, delim_whitespace=True, comment="#").dropna()
    if c95 == True:
        creds = create_95perc_credset(creds, prob_col)
    creds[["chr", "pos", "a1", "a2"]] = creds[SNP_col].str.split(r"[:_]", expand=True)
    creds[["chr", "pos"]] = creds[["chr", "pos"]].astype(np.int32)
    if build == "GRCH38":
        creds = liftover_df(lo, creds, ["chr", "pos"], ["pos_37"])
    if creds is None:
        print(f"\nNo variants in locus {locno} as read from {path_to_credset}")
        return
    GenomicRiskLoci = check_risk_loci(GenomicRiskLoci, creds, locno, filter)
    ref_genes = os.path.join(Ann_path, "ENSG/ENSG.v102.genes.txt")
    if os.path.exists(genes):
        genes = get_genes_in_locus(locno, ref_genes, genes)
        if len(genes) == 0:
            print(f"\nNo genes in locus {locno}")
            return
    else:
        genes = Genes_in_locus_bp_based(locno, GenomicRiskLoci, ref_genes, lo)
    if len(genes) == 0:
        print(f"\nNo genes in locus {locno}")
        return
    if trainset:
        causal = check_TP_genes(genes, trainset)
        if sum(list(genes["ensg"].isin(causal["ensg"]).astype(int))) == 0:
            print(f"\nNo causal genes in locus {outfile}")
            return
    fname = os.path.basename(outfile)
    genes = annotate_credset(
        genes,
        creds,
        build,
        Ann_path,
        cmd_vep,
        vep_cache,
        tabix,
        CADD_file,
        outdir,
        prob_col,
        fname,
        magma=magma_scores,
        trainset=False,
        filter=filter,
    )
    if type(genes) == str:
        print(f"\n{outfile} was not created due to VEP server timeout. Please try again at a later time or set up command line vep")
        return
    genes = POPS_MAGMA_annotation(genes, magma_z, PoPS)
    if trainset:
        causal = check_TP_genes(genes, trainset)
        genes["TP"] = genes["ensg"].isin(causal["ensg"]).astype(int)
    genes = genes.replace(np.nan, 0)
    genes = create_relative_annotation(genes)
    genes = create_locus_metadata(genes, creds, prob_col)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    genes = genes.drop_duplicates(subset=["ensg", "symbol"], keep="first")
    genes.to_csv(
        outfile,
        sep="\t",
        index=False,
    )
    return


# Check if the distance from TP to credset centroid is within 500kb
def check_dist_to_causal(genes, trainset, dist=750000):
    causal = check_TP_genes(genes, trainset)
    if causal["dist"].min() > dist:
        print("\nWARNING: Credible set is too far away from causal gene")
        return False
    return True

def progress_bar(iterable, length=20):
    total = len(iterable)
    progress = 0
    start_time = time.time()

    for i, item in enumerate(iterable):
        yield item
        progress = i + 1
        percentage = progress / total
        progress_length = int(length * percentage)
        elapsed_time = time.time() - start_time

        sys.stdout.write('\033[K')  # Clear the line
        sys.stdout.write(f"\r[{'=' * progress_length}{' ' * (length - progress_length)}] {progress}/{total} ({percentage:.1%}) - Locus processed in: {elapsed_time:.2f}s")
        sys.stdout.flush()
    sys.stdout.write('\n')

# Main fuction to call that will run the annotation for each credset
@runtime
def main(
    credset_file,
    index_file,
    Annotation_dir,
    build,
    pops_out,
    magma_z,
    magma_tissue,
    outdir,
    SNP_col,
    prob_col,
    loci=False,
    genes=False,
    TP=False,
    filter=750000,
    cmd_vep=False,
    vep_cache=False,
    tabix=False,
    CADD_file=False,
    c95=True
):
    # Path for annotation files
    Ann_path = os.path.normpath(Annotation_dir)

    # Load PoPS & MAGMA scores
    PoPS = pd.read_csv(pops_out, sep="\t", engine="pyarrow")
    magma_z = pd.read_csv(magma_z, delim_whitespace=True)
    magma_tissues = pd.read_csv(
        magma_tissue,
        delim_whitespace=True,
        comment="#",
    )
    magma_scores = create_relevance_dict(magma_tissues)

    # Check build
    build = validate_build(build)
    lo = False
    if build == "GRCH38":
        lo = LiftOver("hg38", "hg19")
    GRL_check = "absent"
    # Load GenomicRiskLoci
    if loci != False and loci != None:
        if os.path.exists(loci):
            GenomicRiskLoci = pd.read_csv(
                os.path.normpath(loci), sep="\t", engine="pyarrow"
            )
            GRL_check = "file"
        else:
            GRL_check = "coordinates"
            GenomicRiskLoci = loci

    # Load credset file locations and create arguments
    infiles = []
    if index_file != None:
        credsets = pd.read_csv(index_file, sep="\t")
        if not "GenomicLocus" in credsets.columns:
            credsets["GenomicLocus"] = range(1, len(credsets) + 1)
        if not "Annotfiles" in credsets.columns:
            credsets["Annotfiles"] = None
        for index, row in credsets.iterrows():
            if GRL_check == "absent":
                GenomicRiskLoci = False
                locno = row["GenomicLocus"]
            infiles.append(
                (
                    (row["Filename"], row["GenomicLocus"],row["Annotfiles"]),
                    build,
                    GenomicRiskLoci,
                    Ann_path,
                    outdir,
                    lo,
                    magma_z,
                    PoPS,
                    genes,
                    SNP_col,
                    prob_col,
                    magma_scores,
                    TP,
                    filter,
                    cmd_vep,
                    vep_cache,
                    tabix,
                    CADD_file,
                    c95
                )
            )
    else:
        if GRL_check == "coordinates":
            coordinates = GenomicRiskLoci.split(":")
            GenomicRiskLoci = pd.DataFrame()
            GenomicRiskLoci["GenomicLocus"] = [1]
            GenomicRiskLoci["chr"] = int(coordinates[0].lower().replace("chr", ""))
            GenomicRiskLoci["start"] = int(coordinates[1])
            GenomicRiskLoci["end"] = int(coordinates[2])
            locno = 1
        elif GRL_check == "absent":
            GenomicRiskLoci = False
            locno = 1
        infiles.append(
            (
                (credset_file, locno, None),
                build,
                GenomicRiskLoci,
                Ann_path,
                outdir,
                lo,
                magma_z,
                PoPS,
                genes,
                SNP_col,
                prob_col,
                magma_scores,
                TP,
                filter,
                cmd_vep,
                vep_cache,
                tabix,
                CADD_file,
                c95,
            )
        )
    print('Starting_annotation:')
    for infile in progress_bar(infiles):
        full_annotation_of_credset(*infile)
    return


if __name__ == "__main__":
    main(argv[1:])
