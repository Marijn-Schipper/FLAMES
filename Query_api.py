#import json
import requests
import sys


def query_VEP(chr, pos, a1, a2, build):
    if build.upper() == 'HG19' or build.upper() == 'GRCH37':
        server = "https://grch37.rest.ensembl.org"
    else:
        server = "https://rest.ensembl.org"
    ext = f"/vep/human/region/{chr}:{pos}:{pos}/{a1}/{a2}"
    headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
    r = requests.get(server+ext, headers=headers)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    return decoded

def query_CADD(chr, pos, a1, a2, build):
    if build.upper() == 'HG19' or build.upper() == 'GRCH37':
        server = "https://cadd.gs.washington.edu/api/v1.0/GRCh37-v1.6/"
    else:
        server = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.6/"
    ext = f"{chr}:{pos}"
    headers = { "Content-Type" : "application/json", "Accept" : "application/json"}
    r = requests.get(server+ext, headers=headers)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    CADD = 0
    for var in decoded:
        if var['Ref'] in [a1, a2] and var['Alt'] in [a1, a2]:
            CADD = float(var['RawScore'])

    return CADD