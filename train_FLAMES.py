import pandas as pd
import os
import random
import numpy as np
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GroupShuffleSplit
import pickle as pk
import sys
import annotate as ann
import xgboost as xgb


def load_data(file):
    try:
        df = pd.read_csv(file, sep="\t")
    except:
        print(f"error in loading file{file}")
        sys.exit(f"error in loading file{file}")
    return df


def filter_TP(df, filter=750000):
    dftp = df[df["TP"] == 1.0]
    if int(filter) == 0:
        return df
    if len(dftp) != 1:
        return
    if dftp["distance"].iloc[0] > filter:
        return
    return df


def check_if_duplicate_locus(df, loci_gene_distances):
    ensg = df[df["TP"] == 1.0]["ensg"].iloc[0]
    dist = df[df["TP"] == 1.0]["distance"].iloc[0]
    x = (ensg, dist)
    if x in loci_gene_distances:
        return (loci_gene_distances, False)
    loci_gene_distances.append(x)
    return (loci_gene_distances, True)


def extract_ann_files(annotation_files):
    with open(annotation_files) as f:
        inputfiles = [l.strip() for l in f]
    return inputfiles


def adjust_PoPS_scores(df):
    maxval = max(df["PoPS_Score"].astype(float))
    minval = min(df["PoPS_Score"].astype(float))
    if maxval == minval and maxval > 0.0:
        df["rel_" + "PoPS_Score"] = 1
    elif max(df["PoPS_Score"]) == 0.0:
        df[f"rel_{'PoPS_Score'}"] = df["PoPS_Score"].copy()
    else:
        df[f"rel_{'PoPS_Score'}"] = df["PoPS_Score"].astype(float) - minval
        df[f"rel_{'PoPS_Score'}"] = df[f"rel_{'PoPS_Score'}"].astype(float) / float(
            maxval - minval
        )
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
    ]
    best_predictors = [
        "weighted_distance",
        "TSS_distance",
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
        df = df[df["distance"] <= int(filter)]
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


# Function to train XGB model
def train_XGB(X_train_temp, y_train_temp, seed, group_kfold, groups):
    param_grid = {
        "max_depth": [4, 5, 6],
        "learning_rate": [0.1, 0.01, 0.001],
        "n_estimators": [50, 100, 150, 200],
        "gamma": [0.2, 0.4],
        "subsample": [0.8, 0.9, 1.0],
        "colsample_bytree": [0.5, 0.75],
    }
    XGB_model = xgb.XGBClassifier(
        random_state=seed,
    )
    random_search = RandomizedSearchCV(
        estimator=XGB_model,
        param_distributions=param_grid,
        n_iter=100,
        cv=group_kfold.split(X_train_temp, y_train_temp, groups),
        random_state=seed,
    )
    random_search.fit(X_train_temp, y_train_temp, groups=groups)

    # Get the best model and its hyperparameters
    best_model = random_search.best_estimator_
    best_params = random_search.best_params_
    XGB_model = best_model
    XGB_model.fit(X_train_temp, y_train_temp)
    return XGB_model


def add_into_single_df(list_of_dfs):
    df = pd.concat(list_of_dfs)
    df = df.reset_index(drop=True)
    df["locus"] = pd.factorize(df["filename"])[0]
    return df


def empty_training_loci(df):
    rel_feats = []
    for col in df.columns:
        if (
            "rel_" in col
            and not "MAGMA" in col
            and not "PoPS" in col
            and not "dist" in col
        ):
            rel_feats.append(col)

    df = df[df["TP"] == 1]
    if len(df[df[rel_feats].sum(axis=1) == 0]) >= 1:
        return True
    return False


def main(
    modeldir,
    annotation_files,
    inclusion,
    filter,
    seed,
):
    if str(seed).lower() == "random":
        seed = random.randint(0, 100000)
    command_to_recreate_model = " ".join(sys.argv) + f" --seed {seed}"
    inputfiles = extract_ann_files(annotation_files)
    loci = []
    loci_gene_distances = []
    for f in inputfiles:
        df = load_data(f)
        df["filename"] = f
        if filter != None:
            df = filter_TP(df, filter)
        else:
            df = filter_TP(df)
        if df is None:
            continue
        loci_gene_distances, unique = check_if_duplicate_locus(df, loci_gene_distances)
        if unique is False:
            continue
        if empty_training_loci(df):
            continue
        df = df.drop_duplicates(keep="first", subset=["ensg", "symbol"])
        loci.append(df)
    random.seed(seed)  # Set the seed
    random.shuffle(loci)
    train = add_into_single_df(loci)
    X_train, y_train, df = create_XY_train(train, filter=filter, inclusion=inclusion)
    rel_features = X_train.columns
    groups = df["filename"].values
    gss = GroupShuffleSplit(n_splits=5, train_size=0.8, random_state=42)
    gss.get_n_splits(X_train, y_train, groups)
    print("Training XGB")
    XGB_model = train_XGB(X_train, y_train, seed, gss, groups)
    if not os.path.exists(modeldir):
        os.makedirs(modeldir)
    filename = os.path.join(modeldir, "FLAMES_XGB_model.sav")
    pk.dump(XGB_model, open(filename, "wb"))
    with open(f"{modeldir}/command.txt", "w") as filehandle:
        filehandle.write(f"{command_to_recreate_model}\n")
    with open(f"{modeldir}/features.txt", "w") as filehandle:
        for feature in rel_features:
            filehandle.write(f"{feature}\n")
    with open(os.path.join(modeldir, "train.txt"), "w") as filehandle:
        for f in set(train["filename"]):
            filehandle.write(f"{f}\n")
    return


if __name__ == "__main__":
    main(sys.argv[1:])
