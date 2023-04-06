import annotate_finemapped_fuma_filt as annotate
import os
import sys

filedir = sys.argv[1]
#check if fuma run correctly

if not os.path.exists(os.path.join(filedir, 'genes.txt')):
    print(f'FUMA has not properly run on this folder: {filedir}')
    exit

finemap_output = os.path.join(filedir, 'finemap_output')
finemap_loci = os.listdir(finemap_output)
for f in finemap_loci:
    if 'log_sss' in f:
        log_path = os.path.join(finemap_output, f)
        posteriors = annotate.parse_finemap_logs(log_path)
        if max(posteriors) == posteriors[0]:
            annotate.main(log_path)