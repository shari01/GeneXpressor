from __future__ import annotations
from pathlib import Path
from typing import Optional
import sys
import os
import tempfile
import importlib.resources as pkgres

# rpy2 hardening + console safety (Windows-friendly)
try:
    import rpy2.robjects as ro
    from rpy2.robjects import StrVector
    from rpy2.robjects.packages import importr
    from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
    from rpy2.rinterface_lib import callbacks as rpy2_callbacks
    rpy2_logger.setLevel("ERROR")
except Exception as e:
    print("ERROR: rpy2 is required. Install with: pip install rpy2", file=sys.stderr)
    raise

if hasattr(sys.stdout, "reconfigure"):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

def _safe_console_write(x):
    try:
        s = x.decode("utf-8", errors="replace") if isinstance(x, (bytes, bytearray)) else str(x)
    except Exception:
        s = "[[unprintable rpy2 message]]\n"
    try:
        sys.stdout.write(s)
    except Exception:
        sys.stdout.buffer.write(s.encode(sys.stdout.encoding or "utf-8", errors="replace"))

rpy2_callbacks.consolewrite_print = _safe_console_write
rpy2_callbacks.consolewrite_warnerror = _safe_console_write

# -------- R package bootstrap --------
def ensure_r_packages(pkgs, install_missing: bool = True):
    utils = importr("utils")
    try:
        utils.chooseCRANmirror(ind=1)
    except Exception:
        pass
    installed = [pkg for pkg in ro.r("rownames(installed.packages())")]
    missing = [p for p in pkgs if p not in installed]
    if missing and install_missing:
        print(f"Installing missing R packages: {', '.join(missing)}")
        utils.install_packages(StrVector(missing))

NEEDED_R_PACKAGES = [
    "DESeq2", "BiocParallel", "dplyr", "ggplot2", "ggrepel", "pheatmap",
    "readr", "tidyr", "tibble", "rmarkdown", "RColorBrewer"
]

# -------- load the bundled R script --------
def _load_r_script_text() -> str:
    # resources/deseq2_runner.R is packaged via pyproject tool.hatch.build.include
    with pkgres.as_file(pkgres.files("genexpressor.resources") / "deseq2_runner.R") as p:
        return Path(p).read_text(encoding="utf-8", errors="replace")

_def_loaded = False
def _ensure_r_loaded() -> None:
    global _def_loaded
    if _def_loaded:
        return
    r_code = _load_r_script_text()
    try:
        ro.r(r_code)
    except Exception:
        # Help the user by writing the R to a temp file for line-numbered diagnostics
        tmp_r = os.path.join(tempfile.gettempdir(), "deseq2_runner_tmp.R")
        Path(tmp_r).write_text(r_code, encoding="utf-8", errors="replace")
        ro.r(f'base::source("{tmp_r.replace("\\\\","/")}", echo=FALSE, keep.source=TRUE, encoding="UTF-8")')
    _def_loaded = True

def run_deseq2(
    parent_dir: str,
    pick: str = "AUTO",
    case_level: str = "Disease",
    control_level: str = "Control",
    alpha: float = 0.05,
    lfc_thr: float = 1.0,
    top_labels: int = 20,
    top_heatmap: int = 50,
    make_report: bool = True,
    debug: bool = True,
    threads: int = 1,
    counts: Optional[str] = None,
    metadata: Optional[str] = None,
    dataset_index: Optional[int] = None,
) -> bool:
    """Python wrapper around the bundled R `run_all()` function."""
    # Ensure required R packages exist (do NOT auto-install by default)
    ensure_r_packages(NEEDED_R_PACKAGES, install_missing=False)
    _ensure_r_loaded()

    run_all = ro.globalenv["run_all"]

    parent_dir = str(Path(parent_dir).resolve())
    r_kwargs = dict(
        PARENT_DIR=parent_dir,
        PICK=pick,
        CASE_LEVEL=case_level,
        CONTROL_LVL=control_level,
        ALPHA=float(alpha),
        LFC_THR=float(lfc_thr),
        TOP_LABELS=int(top_labels),
        TOP_HEATMAP=int(top_heatmap),
        MAKE_REPORT=bool(make_report),
        DEBUG=bool(debug),
        threads=int(threads),
        COUNTS_OVERRIDE=(counts if counts else ro.NULL),
        METADATA_OVERRIDE=(metadata if metadata else ro.NULL),
        DATASET_INDEX_OVERRIDE=(int(dataset_index) if dataset_index is not None else ro.NA_Integer),
    )
    ok = run_all(**r_kwargs)
    # ok is an R logical vector; index [0] to get Python bool
    return bool(ok[0])
