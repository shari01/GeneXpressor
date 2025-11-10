from __future__ import annotations
import argparse, sys
from genexpressor.runner import run_deseq2
from genexpressor import __version__

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="genexpressor",
        description="Run DESeq2 analysis via rpy2 using AUTO/ALL dataset discovery or explicit overrides.",
    )
    p.add_argument("--parent_dir", "-p", default=r"C:/Users/shahr/Downloads/Deseq2-pkg",
                   help="Parent directory containing dataset subfolders or files")
    p.add_argument("--pick", default="AUTO",
                   help='Dataset picker: "AUTO", "ALL", or comma-separated indices like "1,3,5"')
    p.add_argument("--case_level", default="Disease")
    p.add_argument("--control_level", default="Control")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--lfc_thr", type=float, default=1.0)
    p.add_argument("--top_labels", type=int, default=20)
    p.add_argument("--top_heatmap", type=int, default=50)
    p.add_argument("--make_report", type=lambda v: str(v).lower() in ("1","true","yes","y"), default=True)
    p.add_argument("--debug", type=lambda v: str(v).lower() in ("1","true","yes","y"), default=True)
    p.add_argument("--threads", type=int, default=1)
    # explicit overrides (single dataset run)
    p.add_argument("--counts", help="Explicit counts file (csv/tsv/txt/xlsx/xls/rds/feather/parquet)")
    p.add_argument("--metadata", help="Explicit metadata file (must have 'condition' and 'sample')")
    p.add_argument("--dataset_index", type=int, help="Force a single dataset index (equivalent to --pick N)")
    p.add_argument("--version", action="version", version=f"GeneXpressor {__version__}")
    return p

def main(argv: list[str] | None = None) -> None:
    argv = argv or sys.argv[1:]
    args = build_parser().parse_args(argv)

    try:
        ok = run_deseq2(
            parent_dir=args.parent_dir,
            pick=args.pick,
            case_level=args.case_level,
            control_level=args.control_level,
            alpha=args.alpha,
            lfc_thr=args.lfc_thr,
            top_labels=args.top_labels,
            top_heatmap=args.top_heatmap,
            make_report=bool(args.make_report),
            debug=bool(args.debug),
            threads=int(args.threads),
            counts=args.counts,
            metadata=args.metadata,
            dataset_index=args.dataset_index,
        )
    except Exception as e:
        print("\nERROR while running DESeq2 pipeline via rpy2:\n", e, file=sys.stderr)
        sys.exit(2)

    if ok:
        print("\n[OK] Completed at least one dataset successfully.")
        sys.exit(0)
    else:
        print("\n[WARN] Pipeline finished & dataset was successfully processed.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
