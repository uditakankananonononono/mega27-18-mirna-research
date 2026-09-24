"""Map miRDB RefSeq transcript accessions to HGNC symbols via MyGene.info batch API.
Output: data/mirdb/refseq2symbol.json (cached; resumable)."""
import json, os, time, urllib.request, urllib.parse
D = os.path.join(os.path.dirname(__file__), '..', 'data', 'mirdb')
out = os.path.join(D, 'refseq2symbol.json')
m = json.load(open(out)) if os.path.exists(out) else {}
ids = [l.strip() for l in open(os.path.join(D, 'nm.txt')) if l.strip() and l.strip() not in m]
for i in range(0, len(ids), 1000):
    b = ids[i:i + 1000]
    data = urllib.parse.urlencode({'q': ','.join(b), 'scopes': 'refseq.rna', 'fields': 'symbol', 'species': 'human'}).encode()
    for att in range(4):
        try:
            r = json.load(urllib.request.urlopen(urllib.request.Request('https://mygene.info/v3/query', data=data), timeout=60))
            break
        except Exception as e:
            time.sleep(3 * (att + 1))
    else:
        raise SystemExit('mygene failed')
    for h in r:
        if 'symbol' in h and h['query'] not in m:
            m[h['query']] = h['symbol']
    for q in b:
        m.setdefault(q, None)
    json.dump(m, open(out, 'w'))
    print(i + len(b), sum(v is not None for v in m.values()), flush=True)
