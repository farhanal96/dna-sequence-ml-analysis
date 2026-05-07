# DNA Sequence ML Analysis

Classify gene types from real human DNA sequences using BioPython and a Random Forest classifier.

## What it does

1. **Fetches** human nucleotide sequences from NCBI (protein-coding mRNA, rRNA, tRNA)
2. **Extracts** biologically meaningful features per sequence:
   - GC content, AT/GC ratio, CpG ratio
   - Mononucleotide, dinucleotide, and trinucleotide (k-mer) frequencies
   - Homopolymer run lengths per nucleotide
3. **Trains** a Random Forest classifier with 5-fold cross-validation
4. **Visualises** confusion matrix, feature importances, and a correlation heatmap

## Project structure

```
dna-sequence-ml-analysis/
├── dna_ml_analysis.ipynb   # Main analysis notebook (run this)
├── dna_features.py         # Feature extraction module
├── requirements.txt        # Python dependencies
├── data/                   # Saved feature matrix (auto-generated)
└── plots/                  # Exported figures (auto-generated)
```

## Quick start

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Launch the notebook
jupyter notebook dna_ml_analysis.ipynb
```

Run all cells top to bottom. NCBI fetching (~450 sequences) takes roughly 2–3 minutes.

## Features extracted

| Feature group | Count | Description |
|---|---|---|
| Sequence stats | 3 | length, GC content, AT/GC ratio |
| CpG ratio | 1 | Observed/expected CpG (promoter signal) |
| Mononucleotide freq | 4 | A, T, G, C proportions |
| Dinucleotide freq | 16 | All 2-mer proportions |
| Trinucleotide freq | 64 | All 3-mer proportions (codon usage) |
| Homopolymer runs | 4 | Max run length per nucleotide |
| **Total** | **92** | |

## Model

- **Algorithm:** Random Forest (300 trees, balanced class weights)
- **Evaluation:** Stratified 5-fold CV + held-out 20 % test set
- **Metrics:** F1-macro, per-class precision/recall, confusion matrix

## Requirements

- Python 3.8+
- See `requirements.txt` for package versions

## NCBI Note

The notebook uses `Entrez.email` as required by NCBI policy. If you have an NCBI API key, uncomment the `Entrez.api_key` line to increase the rate limit from 3 to 10 requests/second.
