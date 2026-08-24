import os
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--filedir', required=True, help="Path to input directory.")
    args = parser.parse_args()
    return args

def make_indexfile(filedir, infile, locus_n):
    
    try:
        outfile = open(f"{filedir}/locus_{locus_n}.cred1", "w")
        print("\t".join(["index", "cred1", "prob1"]), file=outfile)

        
        count = 1
        with open(os.path.join(filedir, infile), "r") as f:
            for line in f:
                if line.startswith("SNP"):
                    continue
                items = line.rstrip("\n").split("\t")
                print("\t".join([str(count), items[0], items[1]]), file=outfile)
                
                
                count += 1
        outfile.close()
    except:
        sys.exit(6)
        
def main():
    args = parse_args()
    filedir = args.filedir
    
    indexfile = open(f"{filedir}/indexfile.txt", "w")
    print("\t".join(["Filename", "GenomicLocus", "Annotfiles"]), file=indexfile)

    locus_n = 1
    for file in os.listdir(filedir):
        if file.endswith(".susie.finemapped"):
            make_indexfile(filedir, file, locus_n)
            # format(os.path.join(filedir, file), locus_n)
            indexfile_row = [f"{filedir}/locus_{locus_n}.cred1", f"{locus_n}", f"{filedir}/annots/annotated_locus_{locus_n}.txt"]
            print("\t".join(indexfile_row), file=indexfile)
            locus_n += 1
    indexfile.close()
    
if __name__ == "__main__":
    main()