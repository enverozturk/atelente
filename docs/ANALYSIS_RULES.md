# Analysis

This guide describes the Shiny workflow and configurations using `cleavage.analysis_scope = "specificity_v3"`.

## Candidate bonds and peptide forms

Each input row is retained as an observation with its modification state. Mapping uses literal sequences and the supplied FASTA, with one-based inclusive coordinates. Alternative protein assignments and repeated placements remain available.

At a bond `k | k+1`:

- A neo-N form starts at `k+1`.
- A neo-C form ends at `k`.
- A spanning form includes residues on both sides.

Candidate termini are identified from incompatibility with the declared digestion rule or, for an internal N terminus, a declared pre-digestion label. Natural protein ends are recorded separately.

Trypsin specificity uses cleavage after K/R except before P; lysine blocking changes this to R except before P. A placement is **specific** when both ends are compatible, **semi-specific** when one is compatible, and **unspecific** when neither is compatible. Reference ends and a possible initiator-Met boundary at position 2 of an M-starting sequence are accepted boundaries for this classification. Label evidence is recorded separately.

Specific and semi-specific placements supply candidate/spanning evidence in this scope. The initial-M bond is reported descriptively. Possible initiator-Met removal is a positional annotation.

In representative-assignment mode, both forms must share a source master accession. The intersection of accession sets supplies the candidate assignments. For each form and shared accession, the first exact occurrence in a single FASTA record must support the same bond. Duplicate records under one accession remain unresolved. Valid shared accessions are alternative mappings of the same measurements.

The overview groups a selected directed form pair, role and comparison across its mappings; one grouped row can therefore have several protein/coordinate assignments.

## Ratios and condition effects

Both forms must have positive finite quantities in the same sample:

```text
sample log2 ratio = log2(neo quantity) - log2(spanning quantity)
ratio log2FC = mean(sample log2 ratio in treatment)
              - mean(sample log2 ratio in control)
```

The displayed neo and spanning effects use those same common samples, so their difference equals the ratio effect. A standalone peptide summary can include additional measurements of that form.

Missing, zero and invalid quantities have separate recorded states. The analysis uses observed quantities without imputation. Forms observed in one condition receive an availability status. Runs without a comparison retain sample quantities with comparison effects/tests unavailable.

## Representative pair selection

Selection is performed for each bond, terminal role and comparison using the following ordered criteria:

1. **Quantitative eligibility.** Common positive measurements must reach `k` in each condition. The default is the smaller of three and the smaller included condition size, with a minimum of one. One-condition observations are handled separately: when only neo is missing from the other condition, spanning must reach `k` there; absence of both components has its own availability state.
2. **End compatibility.** Prefer known digestion-compatible boundaries at the neo's other end and both spanning ends. Reference ends, possible initiator-Met boundaries, unknown metadata and pre-digestion label evidence are neutral for this preference.
3. **Common coverage.** Maximize the smaller condition count, then the total common count.
4. **Other crossed bonds.** Prefer fewer other accepted candidate bonds crossed by the spanning form on the same FASTA record. Exclude the focal bond from this count and retain unknown counts separately.
5. **Pooled joint intensity.** Prefer the larger mean log2 quantity across both forms and all jointly measured included samples.
6. **Stable observation IDs.** Resolve remaining ties reproducibly.

These criteria determine selection independently of effect size and p-values. Internal missed cleavages are excluded from ranking. If all pairs are ineligible, a coverage/ID fallback retains the applicable low-coverage or unavailable status. Alternative pairs remain separate.

Scripted runs can set `representative_selection_k` to a positive integer. Display coverage filters apply to the existing results and leave representative selection fixed.

## Statistical comparisons

An independent-sample OLS model is fitted to sample log2 neo/spanning ratios with limma empirical-Bayes moderation. Testing requires at least three common positive measurements per condition, usable residual variation and at least two eligible pairs in the moderation pool.

Each directed measurement pair is tested once per comparison and test family. Alternative mappings retain that pair's result. Benjamini–Hochberg adjustment uses the tested pairs within each comparison/family.

| Family | Assignment |
| --- | --- |
| `N_standard` | Neo-N with no declared pre-digestion N-label profile |
| `N_TMT_primary` | Neo-N bearing the declared TMT/TMTpro N-terminal label |
| `N_labelled_primary` | Neo-N bearing another declared N-terminal label |
| `N_acetyl_exploratory` | N-acetyl form in a declared label profile |
| `N_unlabelled_exploratory` | Unlabelled N-terminal state in that profile |
| `N_unresolved_exploratory` | Unresolved N-terminal state in that profile |
| `C_secondary` | Neo-C pairs |

Unavailable modification information stays unresolved. The effect describes relative abundance of the measured forms; component abundances and measurement coverage support its interpretation.

## Report summaries

The default display threshold is `|log2FC| >= 1`. Direction summaries use finite, arithmetically valid alternative pairs passing the selected coverage and effect thresholds:

- One direction present: **Increased** or **Decreased**.
- Both directions, with at least 80% of threshold-passing pairs in one direction: **Predominantly increased/decreased**.
- Both directions without that majority: **Mixed directions**.
- Computable pairs below the effect threshold: **Below threshold**.
- No valid finite pair passing coverage: **No computable pairs**.

These summaries count observed pairs. Below-threshold and noncomputable pairs are counted separately. At threshold zero, directions require strictly positive or negative effects.

Pairs with one component unquantified throughout the comparison are retained in Peptides, All pairs and the audit, outside the overview/default Results. Joint measurements in one condition remain in the overview. Other candidates contains the separately reported candidate observations outside the quantitative overview.

Motif panels use fixed representatives for the selected comparison, role and family, with the current coverage and effect thresholds. Increased and decreased ratios have separate panels. Each directed measurement pair contributes once; ambiguous mapping context is omitted at the affected position. Panels report residue composition as percentages with raw counts available and a shared automatic scale.
