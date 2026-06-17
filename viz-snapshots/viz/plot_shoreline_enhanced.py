#!/usr/bin/env python3
"""
CEM shoreline visualization toolkit (enhanced)

Generates multiple figures from a set of ShorePos_*.dat snapshots:
  1) overlay:      all snapshots overlaid (careful when N is large)
  2) subset:       evenly sampled subset overlaid (cleaner)
  3) band:         mean ± std band + first & last curves
  4) heatmap:      time (rows) × alongshore (cols) image
  5) delta:        change relative to a reference snapshot (default: first)

Usage (from a run folder):
  python3 ../../viz/plot_shoreline.py --input "ShorePos_*.dat" --outdir . --ds 100 --subset 8
"""

import argparse, glob, io, re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# --------------------------- IO helpers ---------------------------

def _numeric_key(p: Path):
    """Sort files like ShorePos_2, ShorePos_10 in numeric order."""
    m = re.search(r'(\d+)', p.stem)
    return int(m.group(1)) if m else p.stem

def load_shoreline(path: Path, cols):
    """Load a shoreline snapshot; tolerate single-row and multi-col files."""
    txt = path.read_text(encoding="utf-8", errors="ignore").replace("%", "").strip()
    txt = re.sub(r"\s+", " ", txt)
    arr = np.loadtxt(io.StringIO(txt))

    if arr.ndim == 1:
        y = arr
    else:
        # If multi-column, pick provided y-column (default second col)
        y = arr[:, cols[1]] if arr.shape[1] > 1 else arr[:, 0]
    return y

def load_all(files, cols):
    Ys = []
    for f in files:
        y = load_shoreline(f, cols)
        Ys.append(y)
    Ys = np.vstack(Ys)  # shape: [time, space]
    return Ys

# --------------------------- Plot helpers ---------------------------

def _x_axis(npts, ds):
    x = np.arange(npts)
    return x if ds is None else x * ds

def _fig(figsize):
    return plt.figure(figsize=figsize) if figsize else plt.figure()

def _save(fname, dpi):
    plt.tight_layout()
    plt.savefig(fname, dpi=dpi)
    print(f"Saved {fname}")

# --------------------------- Main plotting ---------------------------

def plot_overlay(files, Ys, ds, show_legend, alpha, figsize, dpi, outpath):
    _fig(figsize)
    x = _x_axis(Ys.shape[1], ds)
    for i, f in enumerate(files):
        plt.plot(x, Ys[i], alpha=alpha, label=f.stem)
    plt.xlabel("Alongshore distance (m)" if ds is not None else "Alongshore cell index")
    plt.ylabel("Shoreline position y")
    plt.title("CEM Shoreline Snapshots (all)")
    if show_legend:
        plt.legend(ncol=2, fontsize=8, framealpha=0.85)
    _save(outpath, dpi)

def plot_subset(files, Ys, ds, k, show_legend, alpha, figsize, dpi, outpath):
    k = max(1, min(k, len(files)))
    idx = np.linspace(0, len(files) - 1, k).astype(int)
    sel_files = [files[i] for i in idx]
    _fig(figsize)
    x = _x_axis(Ys.shape[1], ds)
    for i in idx:
        plt.plot(x, Ys[i], alpha=alpha, label=files[i].stem)
    plt.xlabel("Alongshore distance (m)" if ds is not None else "Alongshore cell index")
    plt.ylabel("Shoreline position y")
    plt.title(f"CEM Shoreline Snapshots (subset k={k})")
    if show_legend:
        plt.legend(ncol=2, fontsize=8, framealpha=0.85)
    _save(outpath, dpi)

def plot_band(files, Ys, ds, figsize, dpi, outpath):
    x = _x_axis(Ys.shape[1], ds)
    mu = Ys.mean(axis=0)
    sd = Ys.std(axis=0)
    _fig(figsize)
    plt.fill_between(x, mu - sd, mu + sd, alpha=0.25, label="mean ± 1σ")
    plt.plot(x, mu, lw=2, label="mean")
    plt.plot(x, Ys[0],  lw=1.2, label=files[0].stem)
    plt.plot(x, Ys[-1], lw=1.2, label=files[-1].stem)
    plt.xlabel("Alongshore distance (m)" if ds is not None else "Alongshore cell index")
    plt.ylabel("Shoreline position y")
    plt.title("CEM shoreline: mean ± std + first/last")
    plt.legend(framealpha=0.85)
    _save(outpath, dpi)

def plot_heatmap(Ys, ds, figsize, dpi, outpath):
    _fig(figsize)
    im = plt.imshow(Ys, aspect="auto", origin="lower", interpolation="nearest")
    plt.colorbar(im, label="Shoreline y")
    plt.ylabel("Snapshot index (time)")
    plt.xlabel("Alongshore distance (m)" if ds is not None else "Alongshore cell index")
    plt.title("CEM shoreline evolution (heatmap)")
    _save(outpath, dpi)

def plot_delta(files, Ys, ds, ref_index, k, show_legend, alpha, figsize, dpi, outpath):
    ref_index = max(0, min(ref_index, Ys.shape[0]-1))
    base = Ys[ref_index]
    x = _x_axis(Ys.shape[1], ds)
    # choose subset
    k = max(1, min(k, Ys.shape[0]))
    idx = np.linspace(0, Ys.shape[0]-1, k).astype(int)
    _fig(figsize)
    for i in idx:
        plt.plot(x, Ys[i] - base, alpha=alpha, label=files[i].stem)
    plt.axhline(0, ls="--", lw=1)
    plt.xlabel("Alongshore distance (m)" if ds is not None else "Alongshore cell index")
    plt.ylabel("Δy (relative to reference)")
    plt.title(f"CEM shoreline change (ref={files[ref_index].stem})")
    if show_legend:
        plt.legend(ncol=2, fontsize=8, framealpha=0.85)
    _save(outpath, dpi)

# --------------------------- CLI ---------------------------

def main():
    p = argparse.ArgumentParser(description="Plot CEM shoreline snapshots (multiple figure types).")
    p.add_argument("--input", default="ShorePos_*.dat", help="Glob for shoreline files")
    p.add_argument("--cols", type=int, nargs=2, default=[0, 1], help="y-column selection if file has multiple columns")
    p.add_argument("--ds", type=float, default=None, help="Alongshore spacing (m). If not set, x is cell index.")
    p.add_argument("--outdir", default=".", help="Directory to write figures")
    p.add_argument("--dpi", type=int, default=220, help="Figure DPI")
    p.add_argument("--figsize", type=float, nargs=2, default=[9, 4], help="Figure size (inches), e.g., --figsize 9 4")
    p.add_argument("--alpha", type=float, default=0.8, help="Line transparency for overlays")
    p.add_argument("--legend", action="store_true", help="Show legends on overlay/subset/delta")
    # subset & delta options
    p.add_argument("--subset", type=int, default=8, help="How many curves to show in subset / delta plots")
    p.add_argument("--ref-index", type=int, default=0, help="Reference index for Δy plot (0=first, -1=last allowed too)")
    # toggles
    p.add_argument("--no-overlay", action="store_true", help="Skip the full overlay figure")
    p.add_argument("--no-subset", action="store_true", help="Skip the subset overlay figure")
    p.add_argument("--no-band", action="store_true", help="Skip the mean±std band figure")
    p.add_argument("--no-heatmap", action="store_true", help="Skip the heatmap figure")
    p.add_argument("--no-delta", action="store_true", help="Skip the Δy (relative change) figure")

    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files = sorted([Path(f) for f in glob.glob(args.input)], key=_numeric_key)
    if not files:
        raise SystemExit(f"No files match: {args.input}")

    Ys = load_all(files, args.cols)  # [time, space]
    figsize = tuple(args.figsize)

    # 1) full overlay
    if not args.no_overlay:
        plot_overlay(files, Ys, args.ds, args.legend, args.alpha, figsize, args.dpi,
                     outdir / "shoreline_overlay.png")

    # 2) subset overlay
    if not args.no_subset:
        plot_subset(files, Ys, args.ds, args.subset, args.legend, args.alpha, figsize, args.dpi,
                    outdir / "shoreline_subset.png")

    # 3) mean ± std band
    if not args.no_band:
        plot_band(files, Ys, args.ds, figsize, args.dpi, outdir / "shoreline_band.png")

    # 4) heatmap
    if not args.no_heatmap:
        plot_heatmap(Ys, args.ds, figsize, args.dpi, outdir / "shoreline_heatmap.png")

    # 5) delta vs reference
    if not args.no_delta:
        ref = args.ref_index if args.ref_index >= 0 else (Ys.shape[0] + args.ref_index)
        plot_delta(files, Ys, args.ds, ref, args.subset, args.legend, args.alpha, figsize, args.dpi,
                   outdir / "shoreline_delta.png")

if __name__ == "__main__":
    main()
