# import modules
import pandas as pd

pd.options.mode.chained_assignment = None  # default='warn'
import os
import multiprocessing as mp
from pyliftover import LiftOver
import numpy as np
from sys import argv
import sys
import time
import Query_api
import subprocess


### Define functions helper functions


# Function to validate build and return a standardized string
def validate_build(build):
    if build.upper() == "HG19" or build.upper() == "GRCH37":
        build = "GRCH37"
    elif build.upper() == "HG38" or build.upper() == "GRCH38":
        build = "GRCH38"
    else:
        sys.exit("Invalid build. Please use build hg19 or hg38")
    return build


# Function to check the range of a genomic location
def check_range(start, stop, loc):
    return start <= loc <= stop


# Wrapper function that print function runtime
def runtime(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        func(*args, **kwargs)
        end = time.time()
        print(f"Runtime: {end - start}")

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
        print(f"Unable to liftover position: {chrom}:{pos}")
    return newpos


# Fuction to liftover a pandas column of genomic locations to different build
def liftover_df(lo, df, incols, outcol):
    print(df)
    df[outcol] = df.apply(
        lambda row: liftover_row(f"{row[incols[0]]}", int(row[incols[1]]), lo),
        axis=1,
    )
    df = df.dropna(subset=[outcol[0]])
    if len(df) == 0:
        print("No variants in credible set could be lifted over")
        return None
    return df


### Functions for annotation of credible sets


# Get the genes in a locus based on the genes file
def get_genes_in_locus(locus_no, path_to_genes_file, GenomicRiskLoci):
    genes = pd.read_csv(path_to_genes_file, sep="\t")
    genes["GenomicLocus"] = genes["GenomicLocus"].astype(str)
    genes["GenomicLocus"] = genes.GenomicLocus.str.split(":")
    genes = genes.explode("GenomicLocus")
    relevant_genes = pd.DataFrame()
    locusname = str(GenomicRiskLoci.iloc[int(locus_no) - 1][0])
    for locus in locusname.split(":"):
        genes_in_loc = genes[genes["GenomicLocus"] == locus]
        relevant_genes = pd.concat([relevant_genes, genes_in_loc])
    return relevant_genes


# Get the genes in locus based on physical location of the locus
def Genes_in_locus_bp_based(locus_no, GenomicRiskLoci, baseline_genes, lo=None):
    Locus = GenomicRiskLoci[
        GenomicRiskLoci["GenomicLocus"].astype(int) == int(locus_no)
    ]
    print(Locus)
    start = Locus["start"].iloc[0]
    end = Locus["end"].iloc[0]
    chrom = Locus["chr"].iloc[0]

    if chrom == "X":
        chrom = "23"
    if not lo == None:
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
    if len(TP_genes) > 500:
        genes_in_locus["TP"] = genes_in_locus["symbol"].isin(TP_genes)
    else:
        genes_in_locus = genes_in_locus[
            genes_in_locus["ensg"].isin(TP_genes)
            | genes_in_locus["symbol"].isin(TP_genes)
        ]
    return genes_in_locus


# Calculate the distances to the TSS and gene body from the credible set
def get_TSS_distances(creds, genes_in_locus, build):
    distances = []
    weighted_distances = []
    TSS_distances = []
    weighted_TSS_distances = []
    poscol = "pos_37" if "pos_37" in creds.columns else "pos"
    centroid = sum(creds[poscol] * creds["prob1"]) / sum(creds["prob1"])
    most_likely_snp = creds.loc[creds["prob1"].idxmax()][poscol]
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
def get_GTEX_eQTLs(creds, genes, finemapped_gtex, build, magma_scores=False):
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
        genepath = os.path.join(finemapped_gtex, row["ensg"].split(".")[0] + ".txt")
        if not os.path.isfile(genepath):
            CLPP_sum.append(sumval)
            CLPP_max.append(maxval)
            tissue_specific_CLPP_sum.append(sumval)
            tissue_specific_CLPP_max.append(maxval)
            continue
        eqtls = pd.read_csv(os.path.join(genepath), sep="\t", engine="pyarrow")

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
            sumval = sum(list(overlap["Probability"] * overlap["prob1"]))
            maxval = max(list(overlap["Probability"] * overlap["prob1"]))
        CLPP_sum.append(sumval)
        CLPP_max.append(maxval)
        if type(magma_scores) == dict and overlap.shape[0] > 0:
            overlap["TISSUE"] = overlap["TISSUE"].str.slice(stop=27)
            overlap["scores"] = overlap["TISSUE"].map(magma_scores)
            overlap["Probability"] = overlap["Probability"] * overlap["scores"]
            sumval = sum(list(overlap["Probability"] * overlap["prob1"]))
            maxval = max(list(overlap["Probability"] * overlap["prob1"]))
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
def get_eqtlgen_eQTLs(creds, genes, eQTL_dir, build):
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
            f"{eQTL_dir}", f"chr{chrom}", "post", f"{row['ensg'].upper()}_new.out_post"
        )
        if not os.path.isfile(genepath):
            CLPP_sum.append(sumval)
            CLPP_max.append(maxval)
            continue
        eqtls = pd.read_csv(os.path.join(genepath), sep="\t", engine="pyarrow")
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
            sumval = sum(list(overlap["Causal_Post._Prob."] * overlap["prob1"]))
            maxval = max(list(overlap["Causal_Post._Prob."] * overlap["prob1"]))
        CLPP_sum.append(sumval)
        CLPP_max.append(maxval)
    genes["CLPP_eQTLgen_sum"] = CLPP_sum
    genes["CLPP_eQTLgen_max"] = CLPP_max
    return genes


def get_Promoters(prom_dir, genes, creds, build):
    sums = []
    maxs = []
    proms = pd.read_csv(
        os.path.join(prom_dir, f"Proms_per_transcript_{build}.txt"), sep="\t"
    )
    for gene in genes["ensg"]:
        prom_df = proms[proms["ensg"] == gene]
        if len(prom_df) == 0:
            sums.append(0)
            maxs.append(0)
            continue
        ## create a matrix that maps which variant is in which captured region
        matrix = create_positional_mapping_matrix(prom_df, creds)
        matrix2 = pd.DataFrame(np.outer(np.ones(len(prom_df)), creds["prob1"]))
        matrix2 = matrix2.where(matrix, other=0)
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["Promoter_sum"] = sums
    genes["Promoter_max"] = maxs
    return genes


# Get the GTEx eQTLs in a locus
def get_VEP(creds, genes, build):
    scores = {}
    VEP_dict = {"HIGH": 1, "MODERATE": 0.6, "LOW": 0.4, "MODIFIER": 0.1}
    VEPs = creds.apply(
        lambda row: Query_api.query_VEP(
            row["chr"], row["pos"], row["a1"], row["a2"], build
        ),
        axis=1,
    )
    if any(VEPs.apply(lambda x: isinstance(x, str))):
        return "Cancel annotation due to timeout in VEP"
    for i, query in VEPs.iteritems():
        variants = query[0]
        PiP = creds["prob1"][i]
        if not "transcript_consequences" in variants.keys():
            continue
        for variant in variants["transcript_consequences"]:
            if not "gene_id" in variant.keys():
                continue
            consequence = float(PiP) * float(VEP_dict[variant["impact"]])
            gene = variant["gene_id"]
            if gene in scores:
                scores[gene].append(consequence)
            else:
                scores[gene] = [consequence]
    genes[["VEP_max", "VEP_sum"]] = genes["ensg"].apply(
        lambda row: pd.Series([max(scores.get(row, [0])), sum(scores.get(row, [0]))])
    )
    return genes


def run_vep_within_environment(
    input_file, output_file, vep_path, build, environment_path
):
    if build.upper() == "GRCH38" or build.upper() == "HG38":
        build = "GRCh38"
    else:
        build = "GRCh37"

    environment_command = f"apptainer exec {environment_path} "
    vep_command = f" {vep_path} -i {input_file} -o {output_file} --cache /home/schipper/.vep --assembly {build} --offline --force_overwrite"
    subprocess.run(environment_command + vep_command, shell=True, check=True)
    return


def cmd_VEP(creds, genes, build, VEP_path, outdir):
    tmp_file_path = os.path.join("/home/schipper/", "tmp_file.vcf")
    output_file_path = os.path.join("/home/schipper/", "vep_output.vcf")
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    environment_path = "/home/schipper/docker_images/ensembl-vep_latest.sif"
    output_lines = creds.copy()
    print(output_lines.head())
    output_lines = output_lines.apply(
        lambda row: [
            f"{row['chr']} {row['pos']} {row['pos']} {row['a1']}/{row['a2']} +",
            f"{row['chr']} {row['pos']} {row['pos']} {row['a1']}/{row['a2']} -",
        ],
        axis=1,
    )
    print(output_lines)
    # Flatten the list of output lines
    output_lines = [line for variant_lines in output_lines for line in variant_lines]

    # Write the lines to a file
    with open(tmp_file_path, "w") as file:
        file.write("\n".join(output_lines))

    run_vep_within_environment(
        tmp_file_path, output_file_path, VEP_path, build, environment_path
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
    vep_df["VEP_weighted"] = vep_df["VEP"] * vep_df["prob1"]
    vep_df["VEP_sum"] = vep_df.groupby("Gene")["VEP_weighted"].transform("sum")
    vep_df["VEP_max"] = vep_df.groupby("Gene")["VEP_weighted"].transform("max")
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
def get_CADD(creds, genes, build):
    creds["CADD"] = creds.apply(
        lambda row: Query_api.query_CADD(
            row["chr"], row["pos"], row["a1"], row["a2"], build
        ),
        axis=1,
    )
    matrix = create_positional_mapping_matrix(genes, creds)
    creds["CADD_weighted"] = creds["CADD"] * creds["prob1"]
    genes["CADD_sum"] = np.dot(matrix, creds["CADD_weighted"])
    max_values = np.nanmax(np.where(matrix, creds["CADD_weighted"], 0), axis=1)
    genes["CADD_max"] = max_values
    return genes


# Get the Jung HiC interaction scores for each gene in the locus
def get_jung_HiC(jung, genes, creds):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(jung, gene + ".txt")):
            sums.append(0)
            maxs.append(0)
            continue
        jungs = pd.read_csv(
            os.path.join(jung, gene + ".txt"), sep=" ", engine="pyarrow"
        )
        ## split jungs 'raw' on '.' and split into three columns called chr, start and end
        jungs[["tmp", "start", "end"]] = jungs["target"].str.split(".", expand=True)
        jungs["p"] = jungs["dist_pvalue"]
        jungs["start"] = jungs["start"].astype(int)
        jungs["end"] = jungs["end"].astype(int)
        jungs = jungs.drop(jungs[jungs["p"] >= 1000].index)
        ## create a matrix that maps which variant is in which captured region
        matrix = create_positional_mapping_matrix(jungs, creds)
        matrix2 = pd.DataFrame(np.outer(jungs["p"], creds["prob1"]))
        matrix2 = matrix2.where(matrix, other=0)
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["Jung_HiC_sum"] = sums
    genes["Jung_HiC_max"] = maxs
    return genes


# Get the Javierre HiC interaction scores for each gene in the locus
def get_javierre_HiC(javierre, genes, creds):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(javierre, gene + ".txt")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        javierres = pd.read_csv(
            os.path.join(javierre, gene + ".txt"), sep="\t", engine="pyarrow"
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
        matrix2 = pd.DataFrame(np.outer(javierres["sum"], creds["prob1"])).where(
            matrix, other=0
        )
        matrix3 = pd.DataFrame(np.outer(javierres["max"], creds["prob1"])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix3.values.max())
    genes["Javierre_HiC_sum"] = sums
    genes["Javierre_HiC_max"] = maxs
    return genes


# Get the ABC EP interaction scores for each gene in the locus
def get_ABC_EP(ABC_EP, genes, creds):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(ABC_EP, gene + ".txt")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        ABCs = pd.read_csv(
            os.path.join(ABC_EP, gene + ".txt"), sep="\t", engine="pyarrow"
        )
        matrix = create_positional_mapping_matrix(ABCs, creds)
        matrix2 = pd.DataFrame(np.outer(ABCs["ABC.Score"], creds["prob1"])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["ABC_EP_sum"] = sums
    genes["ABC_EP_max"] = maxs
    return genes


# Get the ABC CRISPR interaction scores for each gene in the locus
def get_ABC_CRISPR(ABC_CRISPR, genes, creds):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(ABC_CRISPR, gene + ".txt")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        ABCs = pd.read_csv(
            os.path.join(ABC_CRISPR, gene + ".txt"), sep="\t", engine="pyarrow"
        )
        matrix = create_positional_mapping_matrix(ABCs, creds)
        matrix2 = pd.DataFrame(np.outer(ABCs["ABC_max"], creds["prob1"])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["ABC_CRISPR_sum"] = sums
    genes["ABC_CRISPR_max"] = maxs
    return genes


# Get the FANTHOM5 interaction scores for each gene in the locus
def get_F5(F5_ep, genes, creds):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(F5_ep, gene + ".txt")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        F5s = pd.read_csv(os.path.join(F5_ep, gene + ".txt"), sep=" ", engine="pyarrow")
        F5s[["chr", "start", "end"]] = F5s["enhancer"].str.split(":|-", expand=True)
        F5s["start"] = F5s["start"].astype(int)
        F5s["end"] = F5s["end"].astype(int)
        matrix = create_positional_mapping_matrix(F5s, creds)
        matrix2 = pd.DataFrame(np.outer(F5s["correlation"], creds["prob1"])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["F5_ep_sum"] = sums
    genes["F5_ep_max"] = maxs
    return genes


# Get Genehancer predictions for each gene in locus
def get_Genehancer(genehancer, genes, creds):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(genehancer, gene + ".txt")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        genehancer_scores = pd.read_csv(
            os.path.join(genehancer, gene + ".txt"), sep="\t", engine="pyarrow"
        )
        matrix = create_positional_mapping_matrix(genehancer_scores, creds)
        matrix2 = pd.DataFrame(
            np.outer(genehancer_scores["score"], creds["prob1"])
        ).where(matrix, other=0)
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["GeneHancer_max"] = maxs
    genes["GeneHancer_sum"] = sums
    return genes


# Get Roadmap Epigenomics scores for each gene in locus
def get_RoadmapEpi(Roadmap, genes, creds):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(Roadmap, gene + ".txt")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        Epi = pd.read_csv(
            os.path.join(Roadmap, gene + ".txt"), sep="\t", engine="pyarrow"
        )
        matrix = create_positional_mapping_matrix(Epi, creds)
        matrix2 = pd.DataFrame(np.outer(Epi["max"], creds["prob1"])).where(
            matrix, other=0
        )
        sums.append(matrix2.values.sum())
        maxs.append(matrix2.values.max())
    genes["RoadmapEpi_sum"] = sums
    genes["RoadmapEpi_max"] = maxs
    return genes


# Get EpiMap scores for each gene in locus
def get_EpiMap(EpiMap, genes, creds):
    sums = []
    maxs = []
    for gene in genes["ensg"]:
        if not os.path.exists(os.path.join(EpiMap, gene + ".txt")):
            sums.append(0.0)
            maxs.append(0.0)
            continue
        Epi = pd.read_csv(
            os.path.join(EpiMap, gene + ".txt"), sep="\t", engine="pyarrow"
        )
        matrix = create_positional_mapping_matrix(Epi, creds)
        matrix2 = pd.DataFrame(np.outer(Epi["corr"], creds["prob1"])).where(
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
    outdir,
    magma=False,
    trainset=False,
    filter=False,
):
    genes = get_TSS_distances(creds, genes, build)
    if trainset:
        if not check_dist_to_causal(genes):
            exit("The true positive causal gene is too far from lead SNP")
    start_time = time.time()
    print("starting_annotation")
    if cmd_vep == False:
        genes = get_VEP(creds, genes, build)
    else:
        genes = cmd_VEP(creds, genes, build, cmd_vep, outdir)
    if type(genes) == str:
        print("Error in querying VEP, exiting")
        return genes
    runtime = time.time() - start_time
    print(f"Runtime of VEP: {runtime} seconds")
    start_time = time.time()
    genes = get_CADD(creds, genes, build)
    runtime = time.time() - start_time
    print(f"Runtime of CADD: {runtime} seconds")
    start_time = time.time()
    genes = get_eqtlgen_eQTLs(
        creds, genes, os.path.join(Ann_path, "eQTLGen", "caviar_farhad"), build
    )
    runtime = time.time() - start_time
    print(f"Runtime of eqtlGEN: {runtime} seconds")
    start_time = time.time()
    genes = get_GTEX_eQTLs(
        creds, genes, os.path.join(Ann_path, "GTEx_v8_finemapping_CAVIAR"), build, magma
    )
    runtime = time.time() - start_time
    print(f"Runtime of GTEx: {runtime} seconds")
    start_time = time.time()
    genes = get_jung_HiC(os.path.join(Ann_path, "Jung_PCHiC"), genes, creds)
    runtime = time.time() - start_time
    print(f"Runtime of Jung: {runtime} seconds")
    start_time = time.time()
    genes = get_javierre_HiC(os.path.join(Ann_path, "Javierre_PCHiC"), genes, creds)
    runtime = time.time() - start_time
    print(f"Runtime of Javierre: {runtime} seconds")
    start_time = time.time()
    genes = get_ABC_EP(os.path.join(Ann_path, "ABC_EP"), genes, creds)
    genes = get_ABC_CRISPR(os.path.join(Ann_path, "ABC_CRISPR"), genes, creds)
    runtime = time.time() - start_time
    print(f"Runtime of ABC: {runtime} seconds")
    start_time = time.time()
    genes = get_F5(os.path.join(Ann_path, "EP_correlation_F5"), genes, creds)
    genes = get_Genehancer(os.path.join(Ann_path, "GeneHancer"), genes, creds)
    genes = get_RoadmapEpi(os.path.join(Ann_path, "RoadMapEpi"), genes, creds)
    genes = get_EpiMap(os.path.join(Ann_path, "EpiMap"), genes, creds)
    runtime = time.time() - start_time
    print(f"Runtime of rest: {runtime} seconds")
    start_time = time.time()
    genes = get_Promoters(
        os.path.join(Ann_path, "Promoter_regions"), genes, creds, build
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
        # if any(x in feature for x in ["max", "sum", "PoPS", "MAGMA_Z", "dist"]):
        #     scaler = RobustScaler()
        #     data = genes[feature]
        #     data_normalized = scaler.fit_transform(data.values.reshape(-1, 1))
        #     genes[f"robust_{feature}"] = data_normalized
    return genes


# Process entire locus
def full_annotation_of_credset(
    path_to_credset,
    build,
    GenomicRiskLoci,
    Ann_path,
    outdir,
    lo,
    magma,
    PoPS,
    genes,
    magma_scores=False,
    trainset=False,
    filter=False,
    cmd_vep=False,
):
    # load credset
    locno = path_to_credset[1]
    loc_out = locno
    if len(path_to_credset[0].split(".")) > 2:
        loc_out = f'{locno}.{path_to_credset[0].split(".")[-2]}'
    path_to_credset = path_to_credset[0]
    print(locno)
    print(path_to_credset)
    outfile = os.path.join(outdir, "annotated_credset_" + str(loc_out) + ".txt")
    if os.path.exists(outfile):
        print("Already annotated")
        return
    creds = pd.read_csv(path_to_credset, delim_whitespace=True, comment="#")
    creds = creds.dropna()
    if filter:
        # Sort the dataset by the 'probability' column in descending order
        sorted_data = creds.sort_values("prob1")

        # Calculate the cumulative sum of the probabilities
        sorted_data["cumulative_sum"] = sorted_data["prob1"].cumsum()

        # Filter the DataFrame based on the cumulative sum condition
        creds = sorted_data[sorted_data["cumulative_sum"] >= 0.05]

        # Reset the index of the filtered DataFrame
        creds.reset_index(drop=True, inplace=True)
        creds = creds.drop(columns=["cumulative_sum"])

    creds[["chr", "pos", "alleles"]] = creds["cred1"].str.split(":", expand=True)
    creds[["chr", "pos"]] = creds[["chr", "pos"]].astype(np.int32)
    creds[["a1", "a2"]] = creds["alleles"].str.split("_", expand=True)

    # liftover variants if build is GRCH38

    if build == "GRCH38":
        creds = liftover_df(lo, creds, ["chr", "pos"], ["pos_37"])
        if creds is None:
            print(f"No variants in locus {locno} as read from {path_to_credset}")
            return
    # get genes in locus
    if os.path.exists(genes):
        genes = get_genes_in_locus(locno, genes, GenomicRiskLoci)
    else:
        genes = Genes_in_locus_bp_based(
            locno, GenomicRiskLoci, Ann_path + "/ENSG/ENSG.v102.genes.txt", lo
        )

    if trainset:
        causal = check_TP_genes(genes, trainset)
        if sum(list(genes["ensg"].isin(causal["ensg"]).astype(int))) == 0:
            return
    genes = annotate_credset(
        genes,
        creds,
        build,
        Ann_path,
        cmd_vep,
        outdir,
        magma=magma_scores,
        trainset=False,
        filter=filter,
    )
    if type(genes) == str:
        print(f"{outfile} was not created due to VEP timeout")
        return
    # genes = POPS_MAGMA_annotation(genes, magma, PoPS)
    if trainset:
        causal = check_TP_genes(genes, trainset)
        genes["TP"] = genes["ensg"].isin(causal["ensg"]).astype(int)
    genes = genes.replace(np.nan, 0)
    print(outfile)
    genes = genes[genes["distance"] <= 750000]
    genes = create_relative_annotation(genes)
    genes["genes_in_locus"] = len(genes)
    genes["credset_size"] = len(creds)
    genes["genes_within_50kb"] = len(genes[genes["distance"] <= 50000])
    genes["genes_within_100kb"] = len(genes[genes["distance"] <= 100000])
    genes["genes_within_250kb"] = len(genes[genes["distance"] <= 250000])
    genes["Highest_PIP"] = creds["prob1"].max()
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    genes.to_csv(
        os.path.join(outfile),
        sep="\t",
        index=False,
    )
    return


# Check if the distance from TP to credset centroid is within 500kb
def check_dist_to_causal(genes, trainset, dist=500000):
    causal = check_TP_genes(genes, trainset)
    if causal["dist"].min() > dist:
        print("WARNING: Credible set is too far away from causal gene")
        return False
    return True


# Main fuction to call that will run the annotation for each credset
@runtime
def main(
    credset_dir,
    Annotation_dir,
    build,
    pops_out,
    MAGMA_outdir,
    outdir,
    loci,
    genes,
    TP=False,
    filter=False,
    processors=False,
    cmd_vep=False,
):
    # Path for annotation files
    Ann_path = os.path.normpath(Annotation_dir)

    # Load GenomicRiskLoci
    GenomicRiskLoci = pd.read_csv(os.path.normpath(loci), sep="\t", engine="pyarrow")

    # Load PoPS & MAGMA scores
    PoPS = False
    magma_scores = False
    magma = False
    PoPS = pd.read_csv(pops_out, sep="\t", engine="pyarrow")
    magma = pd.read_csv(MAGMA_outdir + "magma.genes.out", sep="\t", engine="pyarrow")
    magma_tissues = pd.read_csv(
        MAGMA_outdir + "magma_exp_gtex_v8_ts_avg_log2TPM.txt.gsa.out",
        delim_whitespace=True,
        comment="#",
    )
    magma_scores = create_relevance_dict(magma_tissues)

    # Check build
    build = validate_build(build)
    if build == "GRCH38":
        lo = LiftOver("hg38", "hg19")
    else:
        lo = False

    # Load credset file locations and create arguments for multiprocessing
    infiles = []
    if os.path.isdir(credset_dir):
        for file in os.listdir(credset_dir):
            if file.endswith(".cred") or file.endswith(".cred1"):
                splitfile = file.split(".")
                locno = splitfile[0]
                locno = locno.split("_")[1]
                infiles.append(
                    (
                        (os.path.join(credset_dir, file), locno),
                        build,
                        GenomicRiskLoci,
                        Ann_path,
                        outdir,
                        lo,
                        magma,
                        PoPS,
                        genes,
                        magma_scores,
                        TP,
                        filter,
                        cmd_vep,
                    )
                )
    else:
        with open(credset_dir, "r") as f:
            for l in f:
                splitline = l.split("\t")
                infiles.append(
                    (
                        (splitline[0], splitline[1].strip()),
                        build,
                        GenomicRiskLoci,
                        Ann_path,
                        outdir,
                        lo,
                        magma,
                        PoPS,
                        genes,
                        magma_scores,
                        TP,
                        filter,
                        cmd_vep,
                    )
                )
    ###THIS IS FOR LATER IMPLEMENTATION OF MULTIPROCESSING
    # if processors:
    #     pool = mp.Pool(processes=int(processors))
    # else:
    #     pool = mp.Pool()
    # apply the function to each file in parallel, passing additional arguments
    # pool.starmap(full_annotation_of_credset, infiles)

    # # wait for the worker processes to complete
    # pool.close()
    # pool.join()
    for infile in infiles:
        full_annotation_of_credset(*infile)
    return


if __name__ == "__main__":
    main(argv[1:])
