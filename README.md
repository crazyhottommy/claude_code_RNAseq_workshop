# RNAseq with Claude Code — a 1-hour workshop

Analyze a real bulk-RNAseq dataset **end-to-end by prompting [Claude Code](https://claude.com/claude-code)**: from salmon quantifications to differential expression, visualization, and GO enrichment.

**📖 Workshop website:** https://crazyhottommy.github.io/claude_code_RNAseq_workshop/

The dataset is [GSE197576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197576) — human cells, **hypoxia vs normoxia**, 4 single-end samples. The repo ships the real salmon `quant.sf` files so you can start analyzing immediately.

## Quickstart

```bash
git clone https://github.com/crazyhottommy/claude_code_RNAseq_workshop.git
cd claude_code_RNAseq_workshop

Rscript install.R          # R/Bioconductor deps (tximport, DESeq2, clusterProfiler, …)
claude                     # launch Claude Code in the repo and follow the website
```

Reproduce every result in one shot:

```bash
Rscript scripts/analysis.R   # -> results/de_hypoxia_vs_normoxia.csv + results/figures/
```

## What's in here

| Path | What it is |
|------|------------|
| `*.qmd` | The Quarto workshop website (rendered to `docs/`) |
| `data/salmon/*/quant.sf` | Real salmon quantifications, 4 samples |
| `data/samples.csv` | Sample sheet (condition assignments) |
| `data/tx2gene.csv`, `data/gene_name_map.csv` | GENCODE v45 transcript→gene + gene→symbol maps |
| `scripts/01–03_*.sh` | Preprocessing: FASTQ download → salmon index → quant |
| `scripts/make_tx2gene.R` | Build the tx2gene / gene-name maps from the GTF |
| `scripts/analysis.R` | Full reference downstream analysis |
| `CLAUDE.md` | Project context that primes Claude Code |
| `environment.yml`, `install.R` | Conda (salmon) + R/Bioconductor dependencies |

## Regenerate the data from raw FASTQ (optional)

```bash
conda env create -f environment.yml && conda activate rnaseq
MODE=SUBSET bash scripts/01_download_fastq.sh   # fast subset (how this repo's data was made)
bash scripts/02_build_salmon_index.sh           # GENCODE v45 human index (RAM-heavy)
bash scripts/03_salmon_quant.sh                 # -> data/salmon/<sample>/quant.sf
Rscript scripts/make_tx2gene.R reference/gencode.v45.annotation.gtf.gz
```

## Build the website locally

```bash
quarto render          # outputs to docs/
quarto preview         # live preview
```

GitHub Pages serves from the `main` branch `/docs` folder. Executed R outputs are frozen in `_freeze/`, so publishing never needs to re-run the analysis.

---

Based on two blog posts by [Tommy Tang](https://divingintogeneticsandgenomics.com/): [downstream with tximport + DESeq2](https://divingintogeneticsandgenomics.com/post/downstream-of-bulk-rnaseq-read-in-salmon-output-using-tximport-and-then-deseq2/) and [preprocessing with salmon](https://divingintogeneticsandgenomics.com/post/how-to-preprocess-geo-bulk-rnaseq-data-with-salmon/).
