import pandas as pd
from sys import argv
import pickle as pk
import os
import train_FLAMES as tf
import optimize_FLAMES as optf


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
    df = tf.add_into_single_df(dfs)
    # Will be removed later if not used for benchmarking
    if not "TP" in list(df.columns):
        df["TP"] = 1
    X_val, y_val, df = tf.create_XY_train(df, filt)
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
            "raw_FLAMES_score",
            "FLAMES_score",
        ]
    ]
    print(f"output written to {os.path.join(outdir, name)}.txt")
    df.to_csv(f"{outdir}/{name}.txt", sep="\t", index=False)
    return


if __name__ == "__main__":
    main(argv[1:])
