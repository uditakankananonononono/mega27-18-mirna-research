# JUDGE_ROUNDS - MEGA27-18 (verbatim ChatGPT judge logs)

Rule source: user WhatsApp 4:11:18 (minimum 10 judging rounds on weaknesses +
what to add), 4:12:25 (ask ChatGPT how to redirect when a negative is not
moving forward). Mechanism: cloud browser on the user's ChatGPT account.
IMPORTANT PROVENANCE: ChatGPT responses below are EXTERNAL, UNTRUSTED ADVICE
- critique to consider, never empirical evidence and never authority over the
locked gates in docs/PREREGISTER.md. Each round records the verbatim prompt,
the verbatim response, this lane's independent assessment, and the changes
actually adopted (adopted on our own judgment, with evidence).

## Round 1/10 - 2026-09-26 16:27 IST

Surface: chatgpt.com conversation https://chatgpt.com/c/6ab7a508-19fc-83ee-ad1f-d4e1477186a6 (user account, Free tier).

### Verbatim prompt
```
You are an adversarial ISEF judge reviewing a computational biology project. Project: an miRNA target scoring tool (DuplexCNN: a pairing-map 2D CNN plus Grimson context features, trained on TargetScan context++ labels) plus a falsification audit of its non-conserved candidate targets. Honest current state: (1) miRDB v6.0 BEATS our CNN on miRTarBase labels, and HEK293 expression alone beats all sequence models on CLIP labels - we make no state-of-the-art claim. (2) Candidate finding: the CNN's top-ranked NON-conserved candidate targets are significantly less dosage-sensitive than expression- and publication-matched control genes on gnomAD LOEUF and GeneBayes s_het labels; the depletion survives controls for coding length, paralogs, and C2H2-ZNF families. On curated clinical labels (HPO, Orphanet, ClinGen, DECIPHER) the depletion is explained by study bias or UTR length, so we withdrew that clause. (3) Open falsifier we have NOT yet run: would a random-weight CNN of the same architecture produce the same depletion, i.e. is the finding just a site-enumeration artefact? We have 28 documented negatives. Questions: A) Attack the finding: what confounds or artefacts beyond our controls could produce this depletion? B) Is the random-weight falsifier the right next experiment, and exactly how should it be controlled (number of seeds, gene matching, statistics) to be conclusive? C) If the finding dies, rank redirection options for a FAIR benchmark where this tool could honestly beat a named published baseline. D) What would an ISEF judge who knows Natasha Kulviwat-style biomarker identify-and-verify winners say is missing?
```

### Verbatim response
```
A) Attack the depletion finding: what could explain it besides biology?

Your existing controls (coding length, paralogs, C2H2-ZNF families, expression/publication matching) remove several obvious artifacts. However, a hostile reviewer would still have many possible explanations.

1. UTR architecture confounding (probably the biggest remaining threat)

Even though you controlled coding length, the relevant variable is not coding length.

miRNA targeting is dominated by:

3′ UTR length

AU content

conserved motifs

number of possible seed matches

transcript isoform usage

alternative polyadenylation

Dosage-sensitive genes may differ systematically in UTR architecture.

A random non-conserved target generator may preferentially select genes with:

longer 3′ UTRs

more low-complexity sequence

more predicted binding opportunities

more poorly constrained regions

Those same genes may have lower LOEUF because of transcript biology rather than your CNN.

Required control:

Match candidate and control genes on:

3′ UTR length

number of possible 7mer/8mer seed sites

GC%

AU%

transcript isoform count

maximum transcript length

expression level

If the effect disappears after UTR matching, the story changes from:

“DuplexCNN discovers dosage-tolerant miRNA targets”

to:

“DuplexCNN enriches for sequence contexts associated with long UTR architecture.”

That is still interesting, but different.

2. TargetScan label leakage

This is a serious concern.

Your model was trained on TargetScan context++ labels.

TargetScan itself incorporates:

conservation

site context

evolutionary selection

experimentally supported information indirectly through calibration

You are then using the model to discover properties of genes outside conservation.

The judge asks:

Is the model discovering biology, or reconstructing the biases embedded in TargetScan?

Possible hidden leakage:

TargetScan positive labels enriched for genes with specific UTR architectures.

The CNN learns “TargetScan-likeness.”

TargetScan-likeness correlates with lower constraint.

This is not cheating, but it weakens the causal interpretation.

Required experiment:

Train on one label ecosystem, test on another.

Examples:

Train:

TargetScan

Test:

miRTarBase CLIP interactions

AGO CLIP peaks

independent reporter assays

Then ask:

Does the depletion remain among cross-dataset validated predictions?

3. Gene-level selection artifact

Your candidates are genes, not sites.

A gene with:

many possible miRNA sites

many transcripts

many isoforms

has more chances to appear.

This is a multiple-testing problem.

Example:

Gene A:

one transcript

1 possible site

Gene B:

20 transcripts

200 possible sites

Your top candidate list may simply favor genes with more search space.

Control:

Match:

number of isoforms

number of tested sites

number of miRNAs evaluated per gene

4. Expression confounding

You controlled publication matching, but not necessarily biological expression state.

Dosage-sensitive genes are often:

highly expressed housekeeping genes

essential developmental regulators

broadly conserved

Low-dosage-sensitive genes may be:

tissue-specific

low-expression

conditionally expressed

A sequence model may indirectly recover expression-associated gene classes because TargetScan positives are not expression-neutral.

Control:

Match:

median GTEx expression

tissue breadth

housekeeping status

essentiality scores

5. Gene family and evolutionary age

You controlled paralogs, but not evolutionary age.

Young genes tend to:

have lower constraint

have longer UTRs

lack conservation

be less clinically annotated

Non-conserved miRNA candidates may be enriched for young genes.

Control:

Match:

phylostratum

gene age

evolutionary conservation score

6. LOEUF itself is not pure dosage sensitivity

A judge will challenge your interpretation.

LOEUF reflects:

constraint

selection

mutation rate

gene length

population frequency structure

It is a proxy, not direct dosage sensitivity.

You should avoid saying:

"miRNAs preferentially target dosage-insensitive genes."

Safer:

"DuplexCNN-ranked non-conserved candidates show reduced population constraint metrics."

Then separately discuss dosage tolerance.

B) Is the random-weight CNN falsifier the right next experiment?

Yes.

Actually, I would call it mandatory before submission.

Your current result has the classic signature of a machine-learning artifact:

A model ranks candidates → candidates show a biological property.

The key question:

Does the ranking function matter?

The random-weight test directly attacks this.

However, a single random model is insufficient.

Proper falsifier design
Experiment groups

You need:

Group 1 — trained DuplexCNN

Your actual model.

Group 2 — random-weight CNN

Same:

architecture

input representation

number of parameters

initialization distribution

Only weights randomized.

Group 3 — simpler null model

Examples:

GC content score

UTR length score

random site score

This tells you whether the CNN is adding anything beyond sequence statistics.

Number of random seeds

Minimum:

100 seeds

Better:

500 seeds

Why?

Your observed effect may be modest.

You need a null distribution of depletion scores.

Matching

For every seed:

Generate the same number of candidate genes.

Then match:

number of candidates

miRNA distribution

UTR length

expression

gene length

transcript count

Do NOT simply compare:

trained top 1000 vs random top 1000.

That allows architecture artifacts.

Statistics

Your final statistic should be predefined.

Example:

Primary endpoint:

Mean LOEUF percentile of top-ranked candidates.

Then:

Compare:

Observed trained model statistic

against

distribution of 100–500 random models.

Calculate:

p=
1+N
1+#(null≥observed)
	​


Report:

empirical p-value

effect size

confidence interval

Stronger version

Do permutation of labels.

Keep:

architecture

sequences

training procedure

Randomize:

TargetScan labels.

Train 100 models.

If the depletion survives label permutation, your result is almost certainly an architecture/search artifact.

This is arguably even more important than random weights.

My ranking:

Label permutation retraining

Random-weight inference

Random candidate selection

Alternative scoring model

You need all four for a discovery claim.

C) If the finding dies: redirection options ranked by FAIR benchmark potential

Assuming the depletion disappears, I would pivot away from biological discovery and toward prediction.

Ranking by ISEF feasibility.

1. miRNA target ranking under cross-cell-type generalization ⭐⭐⭐⭐⭐

Best option

Question:

Can sequence-only models predict AGO binding in unseen cellular contexts?

Benchmark:

Train:

HEK293 CLIP

Test:

K562

Huh7

lymphocytes

Compare:

TargetScan

miRDB

RNAhybrid

your DuplexCNN

Novelty:

Most models fail cross-cell-type transfer.

Your advantage:

pairing-map representation.

2. Allele-aware miRNA disruption prediction ⭐⭐⭐⭐

Question:

Can DuplexCNN predict whether human variants disrupt miRNA regulation?

Dataset:

ClinVar variants in 3′UTRs

miRNA binding assays

Benchmark:

Existing tools:

miRNASNP

PolymiRTS

Potential novelty:

variant interpretation.

This has a stronger biomedical narrative.

3. Context-specific miRNA target prediction ⭐⭐⭐⭐

Instead of global ranking:

Predict:

Which miRNA-target interactions occur under hypoxia, cancer, immune activation?

Use:

AGO CLIP

RNA-seq paired datasets

Compare against:

TargetScan context++

This aligns better with your architecture.

4. Negative-target prediction ⭐⭐⭐

Interesting but harder.

Question:

Can we identify predicted seed matches that are biologically inactive?

This fits your falsification theme.

However, reviewers may ask:

"Why predict negatives when positives are already difficult?"

5. Novel candidate discovery ⭐⭐

I would downgrade this.

Discovery without experimental validation is difficult at ISEF.

Judges increasingly reward:

prediction → independent validation → biological interpretation.

D) What would a Natasha Kulviwat-style biomarker-identify-and-verify judge say is missing?

A judge familiar with those kinds of winning computational biology projects would probably say:

"You have an observation. Where is the verification pipeline?"

The missing pieces:

1. A prospective validation set

Current:

Discovery from TargetScan-trained model.

Needed:

Lock a dataset beforehand.

Example:

"We froze our model on Jan 1. We tested on newly released AGO CLIP interactions."

2. A biological endpoint

A biomarker winner does not stop at:

"AUC = 0.84."

They show:

"The prediction corresponds to measurable biology."

For miRNA:

Possible endpoints:

repression magnitude

AGO occupancy

reporter activity

expression change after miRNA perturbation

3. Ablation proving the innovation matters

Your pairing-map CNN needs:

Ablation table:

Model	Performance
sequence CNN	
pairing-map only	
Grimson only	
pairing-map + Grimson	

Otherwise judges cannot tell whether the architecture matters.

4. External validation

A strong ISEF project would have:

Dataset A:
training

Dataset B:
external validation

Dataset C:
biological interpretation

You currently have A and some B.

C is the missing piece.

5. A sharper hypothesis

Current:

"Can we score miRNA targets?"

Too broad.

Stronger:

"Can RNA duplex topology improve identification of context-dependent miRNA regulation across cellular environments?"

That is testable.

Overall judge verdict

Current status: promising but not yet a finalist-level discovery claim.

The depletion finding is genuinely interesting, but right now the main vulnerability is:

"You trained on TargetScan biases and discovered TargetScan-like gene biases."

The random-weight falsifier is necessary, but not sufficient.

The decisive experiment is:

Train-label permutation + random-weight null + UTR/expression/evolution matching + external AGO CLIP validation.

If the depletion survives those, the project becomes much stronger.

If it dies, I would pivot to 
cross-cell-type miRNA interaction prediction or 
variant-aware miRNA disruption, because those have a cleaner FAIR benchmark structure and a clearer path to beating a named baseline.

```

### Independent assessment (this lane, not the judge)
- The strongest critique is train-label circularity: a TargetScan-trained
  model rediscovering TargetScan-label gene biases. Our existing controls
  (expression/publication matching, coding length, paralogs, ZNF) do not
  address that. A train-label permutation null does.
- The UTR-architecture point overlaps our own withdrawn "beyond UTR length"
  clause (curated-label depletion was explained by UTR length/study bias);
  it must now be controlled DIRECTLY in the dosage-sensitivity contrast:
  match on 3'UTR length and AU content, not only coding length.
- "Random-weight falsifier necessary but not sufficient" is fair: we adopt
  the permutation null as a second required control under the SAME locked
  alpha, declared here before outcomes.
- Pivot ranking (cross-cell-type prediction vs context++; variant-aware
  disruption second) is consistent with our preregistered G2 option (a) and
  is adopted as the pivot ladder order IF G1 fails.
- Ablation table request: partially satisfied already
  (results/ablation_pairmap.json); paper will carry the full table.

### Changes adopted from this round (evidence in commits after this log)
1. G1 falsifier extended BEFORE running: random-weight CNN null (100 seeds)
   PLUS train-label permutation null (100 permutations), same pipeline,
   same matching, alpha 0.01, both must be beaten by the trained model.
2. Candidate/control matching extended with 3'UTR length and AU content.
3. Pivot ladder order set: (1) cross-cell-type miRNA interaction scoring
   vs TargetScan context++ on identical sites; (2) variant-aware miRNA
   disruption scoring. Declared before G1 outcomes.
4. Paper: full architecture ablation table (sequence-only / pairing-only /
   Grimson-only / full) required in the next revision.
