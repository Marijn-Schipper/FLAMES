import pandas as pd
from sys import argv
import pickle as pk
import os
import annotate as ann
import optimize_FLAMES as optf

def add_into_single_df(list_of_dfs):
    df = pd.concat(list_of_dfs)
    df = df.reset_index(drop=True)
    df["locus"] = pd.factorize(df["filename"])[0]
    return df

def create_XY_train(df, filter=None, inclusion="all"):
    train_features = [
        "genes_in_locus",
        "credset_size",
        "genes_within_50kb",
        "genes_within_100kb",
        "genes_within_250kb",
        "Highest_PIP",
        "rel_MAGMA_Z",
        "MAGMA_Z",
    ]
    best_predictors = [
        "weighted_distance",
        "TSS_distance",
        "weighted_TSS_distance",
        "VEP_max",
        "CADD_sum",
        "CLPP_GTEX_tissue_weighted_eQTL_sum",
        "RoadmapEpi_sum",
        "HACER_PRO-seq_GRO-seq_sum",
        "Promoter_sum",
        "F5_ep_sum",
        "Jung_HiC_sum",
        "GeneHancer_sum",
        "EpiMap_sum",
        "mean_CLPP_eQTL_eQTLCatalog_max",
        "ABC_CRISPR_sum",
        "ABC_EP_sum",
        "max_CLPP_rQTL_eQTLCatalog_max",
        "Javierre_HiC_sum",
        "HACER_CAGE_max",
        "CLPP_eQTLgen_sum",
        "cicero_sum",
    ]
    if filter != None:
        df = df[df["weighted_distance"] <= int(filter)]
    else:
        df = df[df["weighted_distance"] <= 750000]
        df = df.groupby("filename").apply(adjust_PoPS_scores)
    if "gene_biotype" in df.columns:
        df = df[df["gene_biotype"] == "protein_coding"]
    if "type" in df.columns:
        df = df[df["type"] == "protein_coding"]
    df = df.groupby("filename").apply(ann.create_relative_annotation)

    for feature in best_predictors:
        if inclusion.lower() == "all":
            train_features.append(feature)
            train_features.append("rel_" + feature)
        if inclusion.lower() == "rel":
            train_features.append("rel_" + feature)
        if inclusion.lower() == "raw":
            train_features.append(feature)
    na_columns = df[train_features].columns[df[train_features].isna().any()].tolist()
    df[na_columns] = df[na_columns].fillna(0)
    train_X = df[train_features]
    train_y = list(df["TP"])
    return train_X, train_y, df

def main(model_dir, input, pops_combi, filt, outdir, name="FLAMES_scores"):
    for file in os.listdir(model_dir):
        if ".sav" in file:
            model = file
    FLAMES = pk.load(open(model_dir + "/" + model, "rb"))
    infiles = []
    if type(input) == list:
        infiles = input
    elif os.path.isdir(input):
            if "annotated_" in file:
                infiles.append(os.path.join(input, file))
            if len(infiles) == 0:
                raise Exception("No annotated credsets starting with 'annotated_' found in directory")
    else:
        raise Exception("Inputfiles not found or not in the right format")
    dfs = [pd.read_csv(f, sep="\t") for f in infiles]
    for i in range(len(dfs)):
        dfs[i]["filename"] = infiles[i]
    df = add_into_single_df(dfs)
    # Will be removed later if not used for benchmarking
    if not "TP" in list(df.columns):
        df["TP"] = 1
    X_val, y_val, df = create_XY_train(df, filt)
    try:
        with open(model_dir + "/features.txt") as f:
            features = [l.strip() for l in f]
    except FileNotFoundError:
        print(
            "features.txt not found in directory containg model, trying to infer from model"
        )
        features = FLAMES._features
    X_val = X_val[features]
    select_features = [
        feature
        for feature in features
        if "rel" in feature
        and not "PoPS" in feature
        and not "MAGMA" in feature
        and not "dist" in feature
    ]
    for feature in X_val.columns:
        X_val[feature] = X_val[feature].astype(float)
    y_pred = FLAMES.predict_proba(X_val)[:, 1]
    df["XGB_score"] = y_pred
    df = df.groupby("locus").apply(optf.Flames_scoring, pops_combi)
    df = df[
        [
            "locus",
            "filename",
            "symbol",
            "ensg",
            "XGB_score",
            "PoPS_Score",
            "rel_XGB_score",
            "rel_PoPS_Score",
            "FLAMES_raw",
            "FLAMES_scaled",
            "estimated_cumulative_precision",
            "FLAMES_highest",
            "FLAMES_causal",

        ]
    ]
    df['estimated_cumulative_precision'] = df['estimated_cumulative_precision'] * df['FLAMES_highest'] * (df['FLAMES_raw']>0.134).astype(int)
    df.to_csv(f"{outdir}/{name}.raw", sep="\t", index=False)
    df=df[df['FLAMES_causal'] == 1]
    df = df.sort_values("FLAMES_scaled", ascending=False)
    df = df[['locus', 'filename', 'symbol', 'ensg', 'FLAMES_scaled', 'FLAMES_raw', 'estimated_cumulative_precision']]
    df.to_csv(f"{outdir}/{name}.pred", sep="\t", index=False)
    print(f"output written to {os.path.join(outdir, name)}.pred and {os.path.join(outdir, name)}.raw")
    return


if __name__ == "__main__":
    main(argv[1:])
