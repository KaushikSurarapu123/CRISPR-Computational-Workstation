# CRISPR Target & Alignment Workstation

An in-silico guide RNA design and genomic safety evaluation tool built in Python.

---

## Overview

Off-target cleavage is one of the primary safety risks in CRISPR-Cas9 genome editing. This application provides an automated pipeline that isolates 20-base-pair target candidates from raw DNA sequences, evaluates cutting efficiency, and runs real-time chromosomal off-target screening using the non-linear MIT Hsu-Zhang position-weighted matrix model.

---

## Key Features

- Live Reference Genome Streaming: Directly fetches and sanitizes FASTA genomic sequences from the NCBI Entrez API (NC_012920.1).
- MIT Hsu-Zhang Off-Target Scoring: Computes empirical cleavage probabilities using the 20-position weighted mismatch matrix, accounting for seed region sensitivities, inter-mismatch distances, and mismatch counts.
- Cutting Efficiency Engine: Evaluates GC content balance (optimal range: 40-60%) and penalizes premature Pol III transcription termination tracts (TTTT).
- Interactive Visual Analytics: Renders dynamic Plotly scatter plots to compare cutting efficiency against off-target risk across all guide candidates.
- Export-Ready Output: Sorts guide RNAs by a composite global score and provides direct CSV exports for lab use.

---

## Tech Stack

| Component | Technology |
| --- | --- |
| Language | Python 3.10+ |
| User Interface | Streamlit |
| Data & Plotting | Pandas, Plotly Express |
| API / Networking | NCBI Entrez Utilities |

---

## Quick Start

1. **Clone the repository:**
   ```bash
   git clone https://github.com/KaushikSurarapu123/CRISPR-Computational-Workstation.git
   cd CRISPR-Computational-Workstation
   ```

2. Install dependencies:
pip install streamlit pandas plotly

3. Launch the workstation:
streamlit run app.py

---

## Primary Citation

Algorithm implementation modeled after the seminal MIT Hsu-Zhang paper:

Hsu, P., Scott, D., Weinstein, J. et al. DNA targeting specificity of RNA-guided Cas9 nucleases. Nature Biotechnology 31, 827–832 (2013).

---

## License

Distributed under the MIT License. See LICENSE for more information.
