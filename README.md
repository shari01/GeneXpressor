<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>⚙️ GeneXpressor 💡</title>
  <style>
    body {
      font-family: Arial, Helvetica, sans-serif;
      line-height: 1.6;
      margin: 40px;
      background-color: #f8f9fa;
      color: #222;
    }
    h1, h2, h3 {
      color: #2c3e50;
    }
    code, pre {
      background-color: #eee;
      border-radius: 6px;
      padding: 4px 6px;
      font-size: 0.95em;
      color: #333;
    }
    pre {
      padding: 10px;
      overflow-x: auto;
    }
    ul {
      margin-left: 20px;
    }
    .highlight {
      background-color: #e3f2fd;
      padding: 10px;
      border-left: 4px solid #2196f3;
      border-radius: 4px;
    }
    footer {
      margin-top: 40px;
      font-size: 0.9em;
      color: #555;
    }
    a {
      color: #007bff;
      text-decoration: none;
    }
    a:hover {
      text-decoration: underline;
    }
  </style>
</head>
<body>
  <header>
    <h1>⚙️ GeneXpressor 💡</h1>
    <p>
      <strong>GeneXpressor</strong> is a streamlined, user-friendly tool designed for bulk RNA-seq differential expression analysis, built on the robust <code>DESeq2</code> framework.
      It intelligently detects your count and metadata files, automates the entire DESeq2 workflow, and generates publication-quality plots along with a reproducible HTML report—all with a single command.
    </p>
    <p>
      Unlike traditional approaches, GeneXpressor allows you to process multiple datasets simultaneously, eliminating the need to select and analyze each dataset individually, making large-scale analyses faster and more efficient.
    </p>
  </header>

  <section>
    <h2>🚀 Key Features</h2>
    <ul>
      <li>🔍 Auto-discovery of count + metadata files</li>
      <li>🧠 End-to-end DESeq2 analysis with rpy2 integration</li>
      <li>📊 Publication-ready visualizations (volcano, MA, heatmap, QC plots)</li>
      <li>🧾 Optional HTML report summarizing results</li>
      <li>⚙️ Cross-platform (Windows, macOS, Linux)</li>
      <li>🧩 Supports both AUTO and manual dataset selection</li>
      <li>💡 Lightweight setup – no manual R scripting required</li>
    </ul>
  </section>

  <section>
    <h2>💾 Installation</h2>
    <pre><code>pip install GeneXpressor
# or
GeneXpressor --help</code></pre>
  </section>

  <section>
    <h2>🧩 Input Directory Structure</h2>
    <p>Your working directory (<code>--parent_dir</code>) should contain:</p>
    <pre><code>MyProject/
├─ counts.csv         # genes × samples, first column = gene IDs
├─ metadata.csv       # sample info; must include the condition column
└─ (optional extra files)</code></pre>

    <h3>Example metadata.csv:</h3>
    <pre><code>sample   condition
S1       Disease
S2       Control
S3       Disease
S4       Control</code></pre>
  </section>

  <section>
    <h2>⚙️ Run DESeq2 Dependencies in R</h2>
    <pre><code>if(!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")

BiocManager::install(c(
  "DESeq2", "apeglm", "ashr", "ggplot2",
  "pheatmap", "rmarkdown", "knitr"
), ask=FALSE)</code></pre>
  </section>

  <section>
    <h2>🧠 CLI Use Case (Windows PowerShell)</h2>
    <pre><code>genexpressor `
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
  --threads 2</code></pre>
  </section>

  <section>
    <h2>📊 Output Overview</h2>
    <p>After a successful run, GeneXpressor will create a new results folder inside your working directory containing:</p>
    <pre><code>C:\Users\<username>\Downloads\Deseq2-pkg\
├─ results/
│  ├─ deseq2_results.csv
│  ├─ deseq2_significant.csv
│  ├─ volcano_plot.png
│  ├─ ma_plot.png
│  ├─ heatmap.png
│  └─ GeneXpressor_Report.html
└─ logs/
   └─ run.log</code></pre>
  </section>

  <section>
    <h2>🧪 Quick Test</h2>
    <pre><code>python -c "import genexpressor; print(genexpressor.__version__)"
genexpressor --help
# Output:
# 0.1.4
# genexpressor CLI OK</code></pre>
  </section>

  <section>
    <h2>🩺 Troubleshooting</h2>
    <table border="1" cellspacing="0" cellpadding="6">
      <tr>
        <th>Issue</th>
        <th>Cause / Fix</th>
      </tr>
      <tr>
        <td><code>ModuleNotFoundError: genexpressor.cli</code></td>
        <td>Reinstall the package using <code>pip install --upgrade GeneXpressor</code></td>
      </tr>
      <tr>
        <td><code>rpy2 or R_HOME errors</code></td>
        <td>Ensure R is installed and added to PATH (<code>R --version</code>)</td>
      </tr>
      <tr>
        <td><code>DESeq2 not found</code></td>
        <td>Install via <code>BiocManager::install("DESeq2")</code></td>
      </tr>
      <tr>
        <td><code>Permission denied</code></td>
        <td>Run from a directory you own</td>
      </tr>
      <tr>
        <td>Missing plots / report</td>
        <td>Check <code>--make_report true</code> and inspect logs</td>
      </tr>
    </table>
  </section>

  <section>
    <h2>📘 Citation</h2>
    <p>
      If you use <strong>GeneXpressor</strong> in your work, please cite:<br />
      <em>Malik, S. (2025). GeneXpressor: Automated DESeq2 runner for bulk RNA-seq via rpy2.</em> PyPI.  
      <a href="https://pypi.org/project/GeneXpressor" target="_blank">https://pypi.org/project/GeneXpressor</a>
    </p>
  </section>

  <section>
    <h2>👤 Author</h2>
    <p>
      <strong>Sheryar Malik</strong><br />
      Bioinformatics Scientist<br />
      📧 <a href="mailto:sheryarmalik1403@gmail.com">sheryarmalik1403@gmail.com</a><br />
      🔗 GitHub: <a href="https://github.com/shari01/GeneXpressor" target="_blank">shari01 / GeneXpressor</a>
    </p>
  </section>

  <footer>
    <h3>📄 License</h3>
    <p>MIT License © 2025 Sheryar Malik</p>
  </footer>
</body>
</html>
