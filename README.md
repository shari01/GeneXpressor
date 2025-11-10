⚙️ GeneXpressor 💡

GeneXpressor is a streamlined, user-friendly tool designed for bulk RNA-seq differential expression analysis, built on the robust DESeq2 framework. 
It intelligently detects your count and metadata files, automates the entire DESeq2 workflow, and generates publication-quality plots along with a reproducible HTML report—all with a single command. 
Unlike traditional approaches, GeneXpressor allows you to process multiple datasets simultaneously, eliminating the need to select and analyze each dataset individually, making large-scale analyses faster and more efficient.​


🚀 Key Features

🔍 Auto-discovery of count + metadata files

🧠 End-to-end DESeq2 analysis with rpy2 integration

📊 Publication-ready visualizations (volcano, MA, heatmap, QC plots)

🧾 Optional HTML report summarizing results

⚙️ Cross-platform (Windows, macOS, Linux)

🧩 Supports both AUTO and manual dataset selection

 Lightweight setup – no manual R scripting required



Run this in your R console:

if(!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")

BiocManager::install(c(
  "DESeq2", "apeglm", "ashr", "ggplot2",
  "pheatmap", "rmarkdown", "knitr"
), ask=FALSE)


💾 Installation
pip install GeneXpressor
# or

To confirm:
GeneXpressor --help


🧩 Input Directory Structure
Your working directory (--parent_dir) should contain:

MyProject/
├─ counts.csv         # genes × samples, first column = gene IDs
├─ metadata.csv       # sample info; must include the condition column
└─ (optional extra files)


Example metadata.csv:

sample	condition
S1	Disease
S2	Control
S3	Disease
S4	Control

⚙️ CLI USE CASE

Example (Windows PowerShell)

genexpressor `
  --parent_dir "C:\Users\<username>\Downloads\Deseq2-pkg" `
  --pick AUTO `
  --case_level Disease `
  --control_level Control `
  --alpha 0.05 `
  --lfc_thr 2.0 `
  --top_labels 20 `
  --top_heatmap 50 `
  --make_report true `
  --debug true `
  --threads 2


📊 Output Overview
After a successful run, GeneXpressor will create a new results folder inside your working directory containing:

Example output structure:

C:\Users\<username>\Downloads\Deseq2-pkg\
├─ results/
│  ├─ deseq2_results.csv
│  ├─ deseq2_significant.csv
│  ├─ volcano_plot.png
│  ├─ ma_plot.png
│  ├─ heatmap.png
│  └─ GeneXpressor_Report.html
└─ logs/
   └─ run.log



🧪 Quick Test
python -c "import genexpressor; print(genexpressor.__version__)"
genexpressor --help
0.1.4
genexpressor CLI OK

Troubleshooting
Issue	Cause / Fix
ModuleNotFoundError: genexpressor.cli	Reinstall the package 
(pip install --upgrade GeneXpressor)

rpy2 or R_HOME errors	Ensure R is installed and on PATH (R --version)
DESeq2 not found	
Install via BiocManager::install("DESeq2")
Permission denied	
Run from a directory you own
Missing plots / report	Check --make_report true and logs


📘Citation
If you use GeneXpressor in your work, please cite:
Malik, S. (2025). GeneXpressor: Automated DESeq2 runner for bulk RNA-seq via rpy2. PyPI. https://pypi.org/project/GeneXpressor

👤 Author
Sheryar Malik
Bioinformatics Scientist
📧 sheryarmalik1403@gmail.com
🔗 GitHub: shari01 : https://github.com/shari01/GeneXpressor

📄 License
MIT License © 2025 Sheryar Malik