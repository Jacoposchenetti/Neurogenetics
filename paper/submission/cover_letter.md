# Cover letter

*(Template — replace the bracketed fields. Written for Imaging Neuroscience; see the
note at the end for Aperture Neuro / PLOS ONE variants.)*

---

[Date]

The Editors
*Imaging Neuroscience*

Dear Editors,

I am submitting for your consideration the manuscript **"A significant enrichment
that is not: spatial nulls, co-expression, and the imaging transcriptomics of EEG
oscillatory-power genetics."**

The paper reports a null result and a methodological demonstration that I believe are
of direct interest to your readership.

I set out to test a straightforward and previously untested hypothesis: that genes
associated with EEG alpha power — from the ENIGMA-EEG oscillatory-power GWAS — are
transcriptomically concentrated in the cortical generators of the alpha rhythm. Using
the conventional analytic route (MAGMA gene-based statistics, Allen Human Brain Atlas
expression via *abagen*, and a spatial-autocorrelation-preserving spin test), the
hypothesis appeared to be supported: enrichment at *p*_spin = 0.022, reproduced by a
second scoring approach, corroborated by an independent surrogate model
(brainsmash, *p* = 0.018), with the pipeline validated by a positive control that
recovered a known cortical expression gradient at *p*_spin = 2 × 10⁻⁴.

Three further analyses, specified as part of the design rather than added afterwards,
dissolve that conclusion. The enrichment is not band-specific — theta, beta and delta
enrich comparably or more strongly. It does not replicate across alpha phenotypes.
And, decisively, when tested against 10,000 random gene sets of matched size rather
than against rotated brain maps, the alpha gene set is unremarkable
(*p*_geneset = 0.33). The apparent effect is a generic property of how co-expressed
gene sets distribute over this cortical territory — a confound that a spatial-only
null is structurally unable to detect. Nothing survived FDR correction across the
full 24-cell analysis grid, and nominal significance did not survive a change of
parcellation.

I think this is worth reporting for two reasons. The substantive null constrains a
plausible hypothesis at the intersection of EEG genetics and imaging transcriptomics,
a pairing that to my knowledge has not previously been attempted. More generally, the
study is a worked case — with the outcome not known in advance, and with sensitivity
established by positive control — in which a spin-significant, surrogate-corroborated,
mechanistically plausible enrichment turns out to be entirely accounted for by gene
co-expression. The theoretical argument is not new (Fulcher et al., 2021), but a
concrete demonstration of what it costs in practice may be more persuasive to
practitioners than the argument alone. The recommendation that follows is simple and
cheap to adopt: report a gene-set null alongside the spatial null.

All analysis code, derived statistics, per-run logs and the manuscript source are
openly available at https://github.com/Jacoposchenetti/Neurogenetics under an MIT
licence. Raw GWAS and atlas data are not redistributed, but access instructions are
documented. Every figure and number in the manuscript can be regenerated from the
repository.

The manuscript is original, is not under consideration elsewhere, and has not been
published previously. I am an independent researcher; this work received no external
funding, and I declare no competing interests. [If applicable: A preprint is
available at bioRxiv doi:XXXXX.]

I would be glad to suggest reviewers with expertise in imaging transcriptomics and
statistical genetics if that would be helpful.

Thank you for your consideration.

Yours sincerely,

Jacopo Schenetti
Independent Researcher, Rome, Italy
jschenetti@gmail.com · ORCID 0009-0000-9108-3806

---

## Adapting for other venues

**Aperture Neuro.** Shorten to ~2 paragraphs and lead with the methodological
contribution — their remit explicitly covers negative results, methods and
reproducibility. Mention that the full analysis is openly reproducible.

**PLOS ONE.** Emphasise that the study is assessed on methodological rigour rather
than novelty of outcome: state plainly that the primary hypothesis was not supported,
that the design included the falsifying tests from the outset, and that sensitivity
was established by positive control. Add a sentence noting the request for
publication-fee assistance (submitted separately).

**Peer Community In / Peer Community Journal.** Reframe as a request for
recommendation of the preprint rather than a journal submission; the emphasis on open
code and full reproducibility is a strength in that setting.
