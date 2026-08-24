import sys, os

infile = sys.argv[1]

outfile = open(sys.argv[2], "w")
print("\t".join(["index", "cred1", "prob1"]), file=outfile)

count = 1

with open(infile, "r") as f:
    for line in f:
        if line.startswith("SNP"):
            continue
        items = line.rstrip("\n").split("\t")
        print("\t".join([str(count), items[0], items[1]]), file=outfile)
        count += 1