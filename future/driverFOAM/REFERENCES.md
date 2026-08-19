# References and Why They Matter

This is an annotated working bibliography. Peer-reviewed publications, preprints,
official project documentation, and venue guidance are labelled separately so
the manuscript does not accidentally present them as equivalent evidence.

## Scientific workflow and research-software precedents

**R1 — Peer-reviewed.** Köster, J. and Rahmann, S. (2012). “Snakemake—a
scalable bioinformatics workflow engine.” *Bioinformatics*, 28(19), 2520–2522.
[doi:10.1093/bioinformatics/bts480](https://doi.org/10.1093/bioinformatics/bts480).

Why it matters: a successful workflow contribution can be demonstrated within a
domain using readable workflow definitions and workstation/cluster execution;
cross-simulator universality is not a prerequisite.

**R2 — Peer-reviewed.** Di Tommaso, P., Chatzou, M., Floden, E. W., Prieto
Barja, P., Palumbo, E., and Notredame, C. (2017). “Nextflow enables
reproducible computational workflows.” *Nature Biotechnology*, 35, 316–319.
[doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820).

Why it matters: portable execution and reproducibility across compute
environments are central evidence for workflow software.

**R3 — Peer-reviewed.** Deelman, E. et al. (2015). “Pegasus, a workflow
management system for science automation.” *Future Generation Computer
Systems*, 46, 17–35.
[doi:10.1016/j.future.2014.10.008](https://doi.org/10.1016/j.future.2014.10.008).

Why it matters: broad generality claims are supported with heterogeneous
scientific applications/platforms, large workflows, failure handling, data
management, monitoring, and provenance.

**R4 — Peer-reviewed comparative study.** Ahmed, A. E. et al. (2021). “Design
considerations for workflow management systems use in production genomics
research and the clinic.” *Scientific Reports*, 11, 21680.
[doi:10.1038/s41598-021-99288-8](https://doi.org/10.1038/s41598-021-99288-8).

Why it matters: supplies a useful evaluation taxonomy—expressiveness,
modularity, scalability, robustness, reproducibility, interoperability, and
ease of development—and combines a real pipeline with platform/scaling tests.

**R5 — Peer-reviewed.** Adorf, C. S., Dodd, P. M., Ramasubramani, V., and
Glotzer, S. C. (2018). “Simple data and workflow management with the signac
framework.” *Computational Materials Science*, 146, 220–229.
[doi:10.1016/j.commatsci.2018.01.035](https://doi.org/10.1016/j.commatsci.2018.01.035).

**R6 — Peer-reviewed proceedings.** Dice, B. D. et al. (2021). “signac: Data
Management and Workflows for Computational Researchers.” *Proceedings of the
20th Python in Science Conference*, 23–32.
[doi:10.25080/majora-1b6fd038-003](https://doi.org/10.25080/majora-1b6fd038-003).

Why R5–R6 matter: they show how a lightweight Python framework can be published
through a coherent data model, automation, provenance/search benefits, and
real computational-research use cases.

**R7 — Official policy.** Association for Computing Machinery. “Artifact
Review and Badging.”
[ACM policy](https://www.acm.org/publications/policies/artifact-review-and-badging-current).

Why it matters: distinguishes repeatability, reproducibility, and replicability;
the manuscript should use these terms precisely.

**R8 — Peer-reviewed conference paper.** Chen, Z. et al. (2025).
“ScienceAgentBench: Toward Rigorous Assessment of Language Agents for
Data-Driven Scientific Discovery.” *ICLR 2025*.
[Conference paper](https://proceedings.iclr.cc/paper_files/paper/2025/hash/f12b4df26344f3be803c06b555252efe-Abstract-Conference.html).

Why it matters: supports task-level assessment before end-to-end automation
claims, using expert-validated tasks, multiple models/scaffolds, repeated
attempts, executable outputs, scientific results, and cost metrics.

## Agent/OpenFOAM precedents

**R9 — Preprint.** Chen, Y., Zhu, X., Zhou, H., and Ren, Z. (2024).
“MetaOpenFOAM: an LLM-based multi-agent framework for CFD.”
[arXiv:2407.21320](https://arxiv.org/abs/2407.21320).

Why it matters: close OpenFOAM-agent baseline using eight CFD tasks, repeated
testing, pass rate, cost, ablation, and randomness sensitivity. Treat as a
preprint unless a peer-reviewed version is cited.

**R10 — Preprint.** Yue, L., Somasekharan, N., Zhang, T., Cao, Y., and Pan, S.
(2025). “Foam-Agent: An End-to-End Composable Multi-Agent Framework for
Automating CFD Simulation in OpenFOAM.”
[arXiv:2509.18178](https://arxiv.org/abs/2509.18178).

Why it matters: closer precedent for composable agent tools, a larger
multi-physics OpenFOAM task corpus, execution success, comparative baselines,
and HPC/external-mesh capabilities. It is a preprint.

**R11 — Peer-reviewed.** Feng, J., Xu, R., and Chu, X. (2026). “OpenFOAMGPT
2.0: End-to-end, trustworthy automation for computational fluid dynamics.”
*International Journal of Heat and Fluid Flow*, 120, 110399.
[doi:10.1016/j.ijheatfluidflow.2026.110399](https://doi.org/10.1016/j.ijheatfluidflow.2026.110399).

Why it matters: current peer-reviewed OpenFOAM-agent precedent reporting more
than 450 executions. Its case families/repeated parametric executions should be
examined carefully before using the execution count as a generality benchmark.

## solids4foam precedents and official evidence

**R12 — Peer-reviewed software paper.** Cardiff, P., Batistić, I., and Tuković,
Ž. (2025). “solids4foam: A toolbox for performing solid mechanics and
fluid-solid interaction simulations in OpenFOAM.” *Journal of Open Source
Software*, 10(108), 7407.
[doi:10.21105/joss.07407](https://doi.org/10.21105/joss.07407).

Why it matters: direct precedent for publishing a mature OpenFOAM research
toolbox and an obvious external-plugin target.

**R13 — Methods/preprint publication.** Cardiff, P., Karač, A., De Jaeger, P.,
Jasak, H., Nagy, J., Ivanković, A., and Tuković, Ž. (2018). “An open-source
finite volume toolbox for solid mechanics and fluid-solid interaction
simulations.” [arXiv:1808.10736](https://arxiv.org/abs/1808.10736).

Why it matters: demonstrates toolbox architecture and representative numerical
problems with comparisons to finite-element solutions.

**R14 — Official documentation.** solids4foam. “Installing solids4foam from
Source — Testing the Installation.”
[Documentation](https://www.solids4foam.com/installation/installFromSource.html).

Why it matters: explicitly distinguishes smoke tests from reference-value
regression tests.

**R15 — Official documentation.** solids4foam. “Tutorials.”
[Tutorial guide](https://www.solids4foam.com/tutorials/).

Why it matters: documents the solids/fluids/FSI taxonomy, required project
dictionaries, and canonical `Allrun` execution pattern.

**R16 — Official documentation.** solids4foam. “About.”
[Project and citation page](https://www.solids4foam.com/about/).

Why it matters: records the project's OpenFOAM compatibility goals,
single-executable design, implementation philosophy, and preferred citations.

## Research-software principles and venues

**R17 — Official journal information.** OpenFOAM Journal. “About the Journal.”
[Journal scope](https://journal.openfoam.com/index.php/ofj/about).

**R18 — Official author guidance.** OpenFOAM Journal. “Author Guidelines.”
[Author guidelines](https://journal.openfoam.com/index.php/ofj/authorGuidelines).

**R19 — Official submission guidance.** Journal of Open Source Software.
“Submitting a paper to JOSS.”
[JOSS guidance](https://joss.readthedocs.io/en/latest/submitting.html).

**R20 — Official journal information.** Elsevier. “SoftwareX.”
[SoftwareX](https://www.sciencedirect.com/journal/softwarex).

**R21 — Official journal information.** Journal of Open Research Software.
[JORS scope](https://openresearchsoftware.metajnl.com/about).

**R22 — Peer-reviewed principles.** Barker, M. et al. (2022). “Introducing the
FAIR Principles for research software.” *Scientific Data*, 9, 622.
[doi:10.1038/s41597-022-01710-x](https://doi.org/10.1038/s41597-022-01710-x).

Why it matters: provides findability, accessibility, interoperability, and
reusability goals adapted to executable/versioned research software.

## Citation-use notes

- Do not cite R9 or R10 as peer-reviewed evidence.
- Verify final bibliographic metadata against the publisher when preparing the
  manuscript.
- Cite an exact driverFOAM software release and data archive separately from the
  article.
- Cite OpenFOAM and cardiacFoam/solids4foam dependencies according to their own
  requested citation guidance.
- Use the local `references.bib` as a starting point, not as an unquestioned
  source of truth.

