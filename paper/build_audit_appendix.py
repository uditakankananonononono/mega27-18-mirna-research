"""Render compact, honest per-miRNA evidence from existing independent benchmark."""
from pathlib import Path
import json
R=Path(__file__).resolve().parent.parent/'results'
B=json.loads((R/'mirdb_benchmark.json').read_text())
O=Path(__file__).resolve().parent/'benchmark_audit.tex'
assert len(B['per_mirna'])==84 and B['mirdb_wins']==74
lines=[r'\section{Independent benchmark audit: every evaluated miRNA}',
 r'The independent miRTarBase 10.0 label benchmark is the clearest competitive boundary for this method. On 7,909 matched rows, miRDB v6.0 has pooled AUROC 0.7867 versus DuplexCNN 0.6847; miRDB wins on 74 of 84 miRNAs. The tables below preserve each per-miRNA sample count, positive count, and pair of AUROCs from \texttt{results/mirdb\_benchmark.json}. These are existing observations, not 84 independent studies or newly run experiments. The smaller miRNAs have unstable estimates. Winning individual rows cannot reverse the pooled or paired overall result.',
 r'\smallskip']
for block in range(0,84,21):
    lines += [r'\begin{table}[htbp]\centering\scriptsize',r'\begin{tabular}{lrrrrr}\hline',
              r'miRNA & Pairs & Positive & CNN AUROC & miRDB AUROC & $\Delta$ (CNN$-$miRDB) \\ \hline']
    for name,z in list(B['per_mirna'].items())[block:block+21]:
        esc=name.replace('_',r'\_')
        a,b=z['cnn_auroc'],z['mirdb_auroc']
        lines.append(f'{esc} & {z["n"]} & {z["n_pos"]} & {a:.3f} & {b:.3f} & {a-b:+.3f} '+r'\\')
    lines += [r'\hline\end{tabular}',
              f'\\caption{{Independent miRTarBase label audit, miRNAs {block+1}--{min(84,block+21)} of 84, in source-record order. Scores are unweighted per-miRNA AUROC and not a pooled estimate.}}\\label{{tab:mirna-audit-{block//21+1}}}',r'\end{table}']
lines += [r'\subsection{Interpretation of the independent ranking}',
 'The source data include 2,113 positive rows out of 7,909, with miRDB scores present for '+
 '29.8\\% of rows. The originally declared test labels, missing-score handling and '+
 'benchmark construction determine the meaning of the comparison; the per-miRNA list '+
 'must not be filtered after seeing scores. The Wilcoxon paired $p=2.45\\times 10^{-12}$ '+
 'supports the direction of miRDB advantage across the tested miRNAs, not biological '+
 'validation of individual target pairs. Later transfer and random-weight controls '+
 'are separate questions; they cannot rehabilitate a false state-of-the-art claim on this task.']
O.write_text('\n'.join(lines)+'\n')
print('wrote',O)
