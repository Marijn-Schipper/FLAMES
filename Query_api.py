# import json
import sys
import json
import time
import requests
from urllib.parse import urlparse, urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError


###Query CADD
def query_CADD(chr, pos, a1, a2, build):
    if build.upper() == "HG19" or build.upper() == "GRCH37":
        server = "https://cadd.gs.washington.edu/api/v1.0/GRCh37-v1.6/"
    else:
        server = "https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.6/"
    ext = f"{chr}:{pos}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    r = requests.get(server + ext, headers=headers)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded = r.json()
    CADD = 0
    for var in decoded:
        if var["Ref"] in [a1, a2] and var["Alt"] in [a1, a2]:
            CADD = float(var["PHRED"])
    return CADD


###Query VEP
class EnsemblRestClient(object):
    def __init__(self, reqs_per_sec=15):
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, silent, hdrs=None, params=None):
        err = 0
        tries = 0
        if hdrs is None:
            hdrs = {}

        if "Content-Type" not in hdrs:
            hdrs["Content-Type"] = "application/json"

        if params:
            endpoint += "?" + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = Request(self.server + endpoint, headers=hdrs)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1
        except HTTPError as e:
            # check if we are being rate limited by the server
            err = e.code
            if e.code == 429:
                if "Retry-After" in e.headers:
                    retry = e.headers["Retry-After"]
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, silent, hdrs, params)
                elif tries < 10:
                    time.sleep(5)
                    self.perform_rest_action(endpoint, silent, hdrs, params)
                    tries += 1
            elif not silent:
                sys.stderr.write(
                    "Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n".format(
                        endpoint, e
                    )
                )
        return (data, err)

    def get_variants(self, chr, pos, a1, a2):
        tries = 0
        variants, e = self.perform_rest_action(
            endpoint=f"/vep/human/region/{chr}:{pos}:{pos}/{a1}", silent=True
        )
        if e == 400:
            variants, e = self.perform_rest_action(
                endpoint=f"/vep/human/region/{chr}:{pos}:{pos}/{a2}", silent=False
            )
        if e == 429 and tries < 5:
            variants, e = self.perform_rest_action(
                endpoint=f"/vep/human/region/{chr}:{pos}:{pos}/{a1}", silent=True
            )
        if e > 200:
            if e == 429:
                sys.stderr.write(f"error in this query due to timeouts: {chr}_{pos}\n")
                return "Cancel annotation due to timeout"
            else:
                return [{"transcript_consequences": [{"impact": "MODIFIER"}]}]
        return variants


def query_VEP(chr, pos, a1, a2, build):
    client = EnsemblRestClient()
    if build == "GRCH37":
        client.server = "http://grch37.rest.ensembl.org"
    else:
        client.server = "http://rest.ensembl.org"
    variants = client.get_variants(chr, pos, a1, a2)
    return variants
