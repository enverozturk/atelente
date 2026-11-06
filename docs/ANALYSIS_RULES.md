# Analysis rules

These definitions describe the current Shiny workflow and scripted runs using `cleavage$analysis_scope = "specificity_v3"`. Low-level mapping and the separate peptide/protein workflow have distinct entry points. A basic mapping configuration does not automatically enable all rules below.

## Observations, placements and candidate bonds

Each input row remains a distinct observation, including its modification information. Matching uses the supplied FASTA and literal sequences. Coordinates are one-based and inclusive. Alternative protein assignments and repeated sequence placements are retained.

For a bond `k | k+1`, a neo-N candidate starts at `k+1`, a neo-C candidate ends at `k`, and a spanning peptide includes residues on both sides. Candidate termini are supported by incompatibility with the declared digestion rule or, for an internal N terminus, the declared pre-digestion label. Natural protein ends are recorded separately.

The supported profile uses trypsin cleavage after K/R except before P. With lysine blocking it uses R except before P. Specificity is assigned per placement: two compatible boundaries are **specific**, one is **semi-specific**, and neither is **unspecific**. Reference ends and possible initiator-Met removal count as accepted boundaries for this classification. Label evidence remains separate. Unspecific placements do not contribute candidate or spanning evidence in this scope. The initial-M bond remains descriptive and is outside internal pair inference.

In representative-assignment mode, the two forms must share a source master accession. The intersection of their accession sets is used, rather than arbitrarily taking the first accession. For each form and shared accession, the first exact occurrence in a single FASTA record is selected. Both occurrences must support the same bond. A later occurrence does not rescue a pair that fails this rule. Duplicate FASTA records under one accession remain unresolved.

Multiple valid shared accessions are alternative attributions of the same measurements. They do not create independent tests. The overview groups the same selected directed form pair, role and comparison across mappings. A grouped row need not identify a unique protein or coordinate.

Possible initiator-Met removal means an N terminus starting at residue 2 of an M-starting reference. It is a position-based label, not confirmation of processing.

## Sample ratios and comparison effects

A sample ratio requires positive finite quantities for both forms in that same sample:

```text
sample log2 ratio = log2(neo quantity) - log2(spanning quantity)
ratio log2FC = mean(sample log2 ratio in treatment)
              - mean(sample log2 ratio in control)
```

No control/treatment replicate pairing is inferred. The neo and spanning component effects displayed with a pair use the same common samples, so their difference equals ratio log2FC. A standalone peptide effect can use that form's additional measurements and is a different summary.

Missing, zero and invalid quantities remain distinct. Values are not imputed. One-condition observations receive an availability status rather than an invented finite FC. Without a comparison, sample measurements can be shown but FC and p-values are unavailable.

## Representative selection

Selection is separate for each bond, terminal role and comparison. Criteria are applied in order, with later criteria breaking ties. Effect size, p-values and desired directions are not selection criteria.

1. **Quantitative eligibility:** common positive measurements must reach selection `k` in both conditions for an eligible two-condition result. Default `k` is the smaller of three and the smaller included condition size, with a minimum of one. One-condition eligibility is handled separately. If only neo is missing from the other condition, spanning must meet `k` there. Both components absent is a different availability state.
2. **End compatibility:** prefer pairs with the neo's other end and both spanning ends known compatible with digestion. The focal neo end is not penalized. Reference ends, possible initiator-Met boundaries and unknown metadata do not establish this preference. Pre-digestion N-label evidence does not establish a digestion-generated end.
3. **Common coverage:** maximize the smaller condition count, then the total common count.
4. **Other crossed candidate bonds:** prefer fewer other distinct accepted bonds crossed by the spanning form on the same FASTA record. The focal bond is excluded from this count. Unknown is not zero.
5. **Pooled joint intensity:** prefer the larger mean log2 quantity across both forms and all jointly measured included samples.
6. **Stable observation IDs:** resolve the remaining tie reproducibly.

Internal missed cleavages do not rank candidates. If no candidate is eligible, the coverage/ID fallback retains its low-coverage or unavailable status. Alternative pairs are not averaged. Display thresholds and coverage filters do not reselect representatives. This is a predefined reporting rule, not a claim of an optimal estimator.

Scripted runs can set the top-level configuration field `representative_selection_k` to a positive integer. Shiny uses the default described above. This analysis setting is distinct from the report's adjustable display coverage filter.

## Optional tests and families

The pair test fits an independent-sample OLS model to sample log2 ratios with limma empirical-Bayes moderation. It requires at least three common positive measurements per condition, usable residual variation and at least two eligible pairs in the moderation pool. Paired/repeated-measure models are not fitted. These minimum counts are not a power guarantee.

Each directed measurement pair is tested once per family and contrast. BH adjustment is applied to the tested pairs separately within that family and contrast. Mapping links and report filters do not add tests or redefine the BH family.

| Test family | Assignment |
| --- | --- |
| `N_standard` | Neo-N when no pre-digestion N-label profile is declared |
| `N_TMT_primary` | Neo-N with the declared TMT/TMTpro N-terminal label |
| `N_labelled_primary` | Neo-N with another declared N-terminal label |
| `N_acetyl_exploratory` | N-acetyl form in a run with a declared label profile |
| `N_unlabelled_exploratory` | Unlabelled N-terminal state in that profile |
| `N_unresolved_exploratory` | Unresolved N-terminal state in that profile |
| `C_secondary` | Neo-C pairs |

Unavailable modification information does not establish an unlabelled state. A test result describes differential ratios of the measured forms. Attribution to a particular processing mechanism or protease requires additional evidence. Condition-dependent missingness and shared denominators can affect interpretation. Ratios are not cleavage percentages.

## Display summaries

The initial ratio threshold is `|log2FC| >= 1` and can be changed. Direction summaries use coverage-passing finite alternative pairs meeting the positive or negative threshold. When only one threshold-passing direction occurs, the label is **Increased** or **Decreased**. If both directions occur, a share of at least 80% of the threshold-passing pairs in one direction is labelled **Predominantly increased/decreased**. Otherwise the label is **Mixed directions**. If computable pairs exist but none passes the effect threshold, the label is **Below threshold**. If no finite, arithmetically valid pair passes the current coverage filter, it is **No computable pairs**. Below-threshold and noncomputable pairs are counted separately. These are descriptive counts, not independent votes. At threshold zero, the sign must be strictly positive or negative.

Pairs with one component unquantified throughout the comparison are omitted from the overview and default Results, without being moved into Other candidates. Their observations and alternative-pair records remain in Peptides, All pairs and the audit. Joint quantification in only one condition is retained. Other candidates holds the separately reported candidate observations outside the quantitative overview.

Motif panels use the fixed representative pairs for the selected comparison, terminal role and test family. Pairs must have a finite, arithmetically valid ratio effect, pass the current coverage filter and meet the current positive or negative ratio threshold. They enter separate **Ratio increased** and **Ratio decreased** panels. At threshold zero, the sign must be strictly positive or negative.

The panels show observed residue percentages, with raw counts available, on a shared automatic scale for the displayed panels. They use no background correction. A directed measurement pair is counted once per contrast, role and family. Ambiguous mapping context is omitted per position. The display describes composition rather than enrichment.
