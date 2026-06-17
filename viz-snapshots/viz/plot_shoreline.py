#!/usr/bin/env python3
import argparse, glob, io, re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def load_shoreline(path, cols):
    txt = Path(path).read_text(encoding="utf-8", errors="ignore")
    txt = txt.replace("%","").strip()
    txt = re.sub(r"\s+", " ", txt)
    arr = np.loadtxt(io.StringIO(txt))
    if arr.ndim == 1:
        y = arr
    else:
        y = arr[:, cols[1]] if arr.shape[1] > 1 else arr[:,0]
    return y

def main():
    p = argparse.ArgumentParser(description="Plot CEM shoreline snapshots.")
    p.add_argument("--input", default="ShorePos_*.dat", help="Glob for shoreline files")
    p.add_argument("--cols", type=int, nargs=2, default=[0,1], help="Zero-based columns for x,y (used if file has multiple columns)")
    p.add_argument("--ds", type=float, default=None, help="Alongshore spacing (m). If not set, x is cell index.")
    p.add_argument("--labels", default=None, help="Comma-separated legend labels (same order as sorted files)")
    p.add_argument("--out", default="shoreline_snapshots.png", help="Output PNG filename")
    args = p.parse_args()

    files = sorted(glob.glob(args.input))
    if not files:
        raise SystemExit(f"No files match: {args.input}")

    labels = args.labels.split(",") if args.labels else None
    if labels and len(labels) != len(files):
        raise SystemExit("Number of labels must match number of files.")

    xs = None
    plt.figure()
    for i, f in enumerate(files):
        y = load_shoreline(f, args.cols)
        if xs is None:
            n = len(y)
            xs = np.arange(n) if args.ds is None else np.arange(n) * args.ds
        lab = labels[i] if labels else Path(f).stem
        plt.plot(xs, y, label=lab)

    plt.xlabel("Alongshore distance (m)" if args.ds is not None else "Alongshore cell index")
    plt.ylabel("Shoreline position y")
    plt.title("CEM Shoreline Snapshots")
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.out, dpi=220)
    print(f"Saved {args.out}")

if __name__ == "__main__":
    main()
