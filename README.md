# Genomic structural variation across clownfish genomes

Code and small reference tables accompanying **Marcionetti et al.**,
*"Genomic structural variation across clownfish genomes"*.

## What's here

- **`Clownfish_GenomeAssembly_And_SV.ipynb`** — the single master notebook. It
  documents the full analysis, from raw PacBio HiFi data to the structural-variant
  (SV) catalogue and downstream analyses, as a sequence of markdown explanations +
  the exact commands / script invocations used, in the order they were run. Every
  numbered section corresponds to a subsection of the manuscript *Methods*. It is
  **documentation of the workflow, not a script to run top-to-bottom** — most steps
  are SLURM array jobs over large genomic files that are not part of this
  repository.
- **`scripts/`** — every helper script referenced from the notebook, grouped by
  pipeline stage (`1_PacBio_Data_and_Preprocessing/` … `9_Long_Inversions_and_Chrom_Trees/`).
- **`metadata/`** — small, hand-curated input tables consumed throughout the
  pipeline (species lists, chromosome ID/orientation maps, host-category
  assignments). This is reference data the scripts read, not pipeline output —
  computed results (VCFs, alignment files, statistics tables, the SV catalogue)
  are not tracked here; see **Data availability** below.

## Repository layout

```
New_Git/
├── README.md
├── .gitignore
├── Clownfish_GenomeAssembly_And_SV.ipynb
├── metadata/
│   ├── Species.txt / Species_OMA.txt / Species_with_outgroups.txt
│   ├── ChromosomesRename_and_RC.txt / ReverseComplement_Species.txt
│   ├── CLA_vs_BIARef.mapids.txt
│   ├── Host_Categories.txt
│   ├── PacBioSpecies.HOG_Tree.treefile / PacBioSpecies.OnlyClownfish.tree
│   └── BIA_ChromLength.txt / CLA_ChromLength.txt
└── scripts/
    ├── 1_PacBio_Data_and_Preprocessing/
    ├── 2_Genome_Assembly/
    ├── 3_Genome_Quality/
    ├── 4_Genome_Annotation/
    ├── 5_Phylogenomics_OMA_AGORA/
    ├── 6_WGA_and_SV/
    ├── 7_SV_Catalogue/
    ├── 8_SV_and_Host_Specialization/
    └── 9_Long_Inversions_and_Chrom_Trees/
```

Folder numbering follows the notebook's section numbering, which follows the
manuscript's Methods order.

## Study species

16 clownfish genomes newly assembled for this study, complemented with public
assemblies and outgroups (see the notebook's species table for platform/role
details): AKA, AKY, ALL, BIA, CRP, EPH, FRE, LAT, LAZ, MCC, OMA, POL, PRC, PRD,
SAN, SEB (newly assembled); BIARef, CLA, OCE, PRCRef (public references); DTR,
ACH (outgroups).

## Software used

Versions as run for the manuscript; see the notebook for the exact command lines.

| Stage | Tool | Version |
|---|---|---|
| QC / k-mer profiling | FastQC | 0.12.1 |
| | MultiQC | 1.25.1 |
| | Jellyfish | — |
| | GenomeScope | 1.0 |
| Assembly | Hifiasm | 0.19.5 |
| | Flye | 2.9.1 |
| | HiCanu | 2.2 |
| | IPA | 1.8.0 |
| | purge_dups | 1.2.5 |
| | ntLink | 1.3.9 |
| | minimap2 | 2.12 (assembly stage) |
| | BWA | 0.7.17 |
| | samtools | 1.19.2 |
| | P_RNA_scaffolder | — |
| | MUMmer4 | 4.0.0 |
| Assembly QC | BUSCO | — |
| | Inspector | 1.3.1 |
| | Clair3 | — |
| | Sniffles2 | 2.5 |
| | VCFtools | 0.1.16 |
| Annotation | EDTA | 2.1.3 |
| | RepeatModeler2 | — |
| | RepeatMasker | 4.0.7 |
| | BRAKER | 3.0.3 |
| | InterProScan | — |
| | DIAMOND (vs Swiss-Prot / nr) | 2.0.15 |
| | Circos | 0.69-8 |
| Phylogenomics / OMA / AGORA | OMA standalone | — |
| | MAFFT | 7.505 |
| | trimAl | — |
| | AMAS | — |
| | IQ-TREE | 2.2.2.7 |
| | iTOL | 7.5.1 |
| | AGORA | — |
| | syntenyPlotteR | 1.0 |
| WGA & SV calling | minimap2 (`-x asm5`) | 2.12 |
| | MUMmer4 (robustness re-run) | 4.0.0 |
| | SyRI | 1.6.3 |
| | tabix / vcf-merge | 1.21 |
| | plotsr | — |
| SV catalogue | GenomicRanges / rtracklayer (Bioconductor) | — |
| Host specialization | vegan (R, Jaccard/PCoA) | — |
| Long inversions & chromosome trees | IQ-TREE | 2.2.2.7 |
| | ape (R, monophyly / Robinson–Foulds) | — |
| | SAMtools / wally (breakpoint genotyping) | 1.19.2 / 0.6.1 |
| | TopGO | 2.62.0 |

Databases: Swiss-Prot, NCBI **nr** (snapshot 2024-02-07), DFAM Vertebrata.

## Third-party scripts (not vendored)

Referenced by name from the notebook/scripts but not copied into this repo —
install from upstream:

| Tool | Source |
|---|---|
| `pafCoordsDotPlotly.R` / `mummerCoordsDotPlotly.R` | dotPlotly — https://github.com/tpoorten/dotPlotly |
| `AMAS.py` | https://github.com/marekborowiec/AMAS |
| `agora-generic.py`, `src/misc.compareGenomes.py` | AGORA — https://github.com/DyogenIBENS/Agora |

## Data availability

Large inputs/outputs (raw reads, assemblies, BAM/VCF files, `.xlsx` tables,
figures/PDFs) are **not** tracked in this repository (see `.gitignore`) and are
archived separately / available on request. Raw PacBio HiFi data: SRA BioProject
*(accession pending)*.

## SV catalogue

`ClownfishSV_database.html` is a self-contained, interactive catalogue of every
SV (built by `scripts/7_SV_Catalogue/generate_sv_database.py`, DataTables-based)
and is tracked in this repository. If downloaded or cloned locally, just open
it in a browser. It can also be viewed directly at:

https://amarcion.github.io/ClownfishSV/ClownfishSV_database.html

It's a large, self-contained file (~25 MB, all SV data embedded), so it may take
a few seconds to load.


## Citation

If you use this code, please cite Marcionetti et al., *"Genomic structural
variation across clownfish genomes"* (2026).
