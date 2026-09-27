"""Fixed-width absolute UTR-length-bin controls, a post-falsifier descriptive audit.

This is NOT a rescue of the failed random-weight gate. Held-out panel and CNN
ranking come from the existing preregistered gnomAD replication; this control
was specified after those outcomes. Each candidate is matched to a non-candidate
gene in the same absolute UTR-length bin (25, 50, or 100 nt). Three seeded draws; sampling
with replacement within bin because sparse bins cannot supply 200 distinct
non-candidates. No nearest-bin fallback: report and drop unsupported candidates,
then compare CNN and controls on the same matchable subset only. One-sided
sign test is descriptive and does not confer a newly preregistered claim.
"""
import json, random, os
import numpy as np
import pandas as pd
from scipy.stats import binomtest

WIDTH = int(os.environ.get("LENGTH_BIN_WIDTH", "50"))
K = 200
R = 3

def main():
    prior = json.load(open('results/gnomad_replication.json'))
    panel = prior['panel']
    lo = pd.read_csv('data/gnomad/lof_v211.txt.bgz', sep='\t', compression='gzip',
                     usecols=['gene', 'oe_lof_upper']).dropna().groupby('gene').oe_lof_upper.min()
    d = pd.read_csv('results/clip_rows_136.csv', usecols=['mirna','gene','min','len']).drop_duplicates(['mirna','gene'])
    d = d[d.gene.isin(lo.index)]
    rows = {}
    for mi in panel:
        s = d[d.mirna == mi]; cand = s.nsmallest(K,'min').copy()
        pool = s[~s.gene.isin(cand.gene)].copy()
        cand['bin'] = cand.len // WIDTH; pool['bin'] = pool.len // WIDTH
        pools = {int(b): list(x.gene) for b,x in pool.groupby('bin')}
        eligible = cand[cand.bin.isin(pools)].copy()
        excluded = cand[~cand.bin.isin(pools)].copy()
        rng = random.Random(f'{mi}-exact{WIDTH}-length-20260928')
        controls=[]; length_abs_diffs=[]; unique_counts=[]
        for _ in range(R):
            gs = [rng.choice(pools[int(b)]) for b in eligible.bin]
            controls.append(float(np.median(lo.reindex(gs))))
            unique_counts.append(len(set(gs)))
            lengths = pool.set_index('gene').len
            length_abs_diffs.append(float(np.median(np.abs(eligible.len.to_numpy() - lengths.reindex(gs).to_numpy()))))
        rows[mi] = {
            'universe_n':int(len(s)), 'candidate_n':int(len(cand)),
            'matchable_n':int(len(eligible)), 'unmatched_n':int(len(excluded)),
            'unmatched_genes':list(excluded.gene),
            'candidate_median_loeuf_matched_subset':float(np.median(lo.reindex(eligible.gene))),
            'control_median_loeuf':controls,
            'control_mean_loeuf':float(np.mean(controls)),
            'median_absolute_length_difference_bp':length_abs_diffs,
            'sampling_with_replacement':True,
            'unique_control_genes_per_draw':unique_counts,
        }
        print(mi, rows[mi]['matchable_n'], rows[mi]['candidate_median_loeuf_matched_subset'], rows[mi]['control_mean_loeuf'], flush=True)
    wins=sum(v['candidate_median_loeuf_matched_subset']>v['control_mean_loeuf'] for v in rows.values())
    losses=sum(v['candidate_median_loeuf_matched_subset']<v['control_mean_loeuf'] for v in rows.values())
    out={'status':'post-hoc descriptive after failed random-weight falsifier',
         'bin_width_bp': WIDTH, 'draws':R, 'panel_source':'results/gnomad_replication.json',
         'method':'sample with replacement from non-candidates within identical absolute 50-nt bin; exclude candidates in bins with zero pool',
         'sign_test_descriptive':{'wins':wins,'losses':losses,'ties':len(rows)-wins-losses,
                                  'p_one_sided':float(binomtest(wins,wins+losses,0.5,alternative='greater').pvalue) if wins+losses else None},
         'matching_coverage':{'eligible':sum(v['matchable_n'] for v in rows.values()),'total':K*len(rows)},
         'per_mirna':rows}
    with open(f'results/exact_length_bin_{WIDTH}bp.json','w') as fh:json.dump(out,fh,indent=1)
    print(json.dumps({'sign_test':out['sign_test_descriptive'],'coverage':out['matching_coverage']}),flush=True)
if __name__=='__main__':main()
