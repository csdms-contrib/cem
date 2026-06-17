#!/usr/bin/env python3
"""
Animate CEM shoreline snapshots.

Examples (run inside a run folder):
  # basic GIF, 6 fps
  python3 ../../viz/animate_shoreline.py --input "ShorePos_*.dat" --gif shoreline.gif --fps 6

  # slower (interval in ms), with 10 in-betweens between each pair
  python3 ../../viz/animate_shoreline.py --input "ShorePos_*.dat" --gif shoreline_smooth.gif --interval 400 --interp 10

  # export both GIF and MP4, set alongshore spacing to meters, and fix y-limits
  python3 ../../viz/animate_shoreline.py --input "ShorePos_*.dat" --gif s.gif --mp4 s.mp4 --fps 5 --ds 100 --ylim 0 20
"""

import argparse, glob, io, re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# --------------------- helpers ---------------------

def _numeric_key(p: Path):
    m = re.search(r'(\d+)', p.stem)
    return int(m.group(1)) if m else p.stem

def load_shoreline(path: Path, ycol: int = 1):
    """Load one shoreline snapshot (robust to stray % and whitespace)."""
    txt = path.read_text(encoding="utf-8", errors="ignore").replace("%", "").strip()
    txt = re.sub(r"\s+", " ", txt)
    arr = np.loadtxt(io.StringIO(txt))
    if arr.ndim == 1:
        y = arr
    else:
        y = arr[:, ycol] if arr.shape[1] > ycol else arr[:, 0]
    return np.asarray(y, dtype=float)

def build_frames(base_frames, interp: int):
    """Insert linear in-betweens between successive 1D arrays."""
    if interp <= 0:
        return list(base_frames)
    out = []
    for i in range(len(base_frames) - 1):
        a, b = base_frames[i], base_frames[i + 1]
        out.append(a)
        for k in range(1, interp + 1):
            w = k / (interp + 1.0)
            out.append((1 - w) * a + w * b)
    out.append(base_frames[-1])
    return out

# --------------------- main ---------------------

def main():
    ap = argparse.ArgumentParser(description="Animate CEM shoreline snapshots.")
    ap.add_argument("--input", default="ShorePos_*.dat", help="Glob for shoreline files")
    ap.add_argument("--ycol", type=int, default=1, help="y-column if file has multiple columns (0-based)")
    ap.add_argument("--ds", type=float, default=None, help="Alongshore spacing (m); if None, x is cell index")
    ap.add_argument("--gif", default="shoreline.gif", help="GIF filename (set to '' to disable)")
    ap.add_argument("--mp4", default="", help="Optional MP4 filename (requires ffmpeg)")
    ap.add_argument("--fps", type=int, default=None, help="Frames per second (alternative to --interval)")
    ap.add_argument("--interval", type=int, default=None, help="Frame interval in ms (alternative to --fps)")
    ap.add_argument("--interp", type=int, default=0, help="Linear in-between frames per pair (>=0)")
    ap.add_argument("--figsize", type=float, nargs=2, default=[9, 4], help="Figure size (inches)")
    ap.add_argument("--ylim", type=float, nargs=2, default=None, help="y-axis limits, e.g. --ylim 0 20")
    ap.add_argument("--title", default="CEM shoreline animation", help="Figure title")
    args = ap.parse_args()

    files = sorted([Path(f) for f in glob.glob(args.input)], key=_numeric_key)
    if not files:
        raise SystemExit(f"No files match: {args.input}")

    # load base snapshots
    base = [load_shoreline(f, args.ycol) for f in files]
    n = base[0].size
    # x axis
    x = np.arange(n) if args.ds is None else np.arange(n) * args.ds

    # stack for auto y-limits
    stack = np.vstack(base)
    ymin, ymax = stack.min(), stack.max()
    pad = 0.05 * (ymax - ymin if ymax > ymin else 1.0)
    if args.ylim is None:
        ylims = (ymin - pad, ymax + pad)
    else:
        ylims = (args.ylim[0], args.ylim[1])

    # interpolate (optional)
    frames = build_frames(base, max(0, args.interp))

    # fps / interval
    if args.fps is not None and args.interval is not None:
        # prefer fps; if both provided, interval will be ignored
        pass
    if args.interval is None:
        fps = args.fps if args.fps is not None else 6
        interval = int(1000 / max(1, fps))
    else:
        interval = max(1, args.interval)

    # figure
    plt.close("all")
    fig, ax = plt.subplots(figsize=tuple(args.figsize))
    (line,) = ax.plot(x, frames[0], lw=1.6)
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(*ylims)
    ax.set_xlabel("Alongshore distance (m)" if args.ds is not None else "Alongshore cell index")
    ax.set_ylabel("Shoreline position y")
    title = ax.set_title(f"{args.title}  |  {files[0].stem}")

    # map each frame to a title (file-based names; interpolated frames show '→')
    names = []
    for i in range(len(base) - 1):
        names.append(files[i].stem)
        names += [f"{files[i].stem}→{files[i+1].stem}"] * max(0, args.interp)
    names.append(files[-1].stem)

    def update(i):
        line.set_ydata(frames[i])
        title.set_text(f"{args.title}  |  {names[i]}")
        return (line, title)

    anim = FuncAnimation(fig, update, frames=len(frames), interval=interval, blit=True)

    # export GIF
    if args.gif:
        writer = PillowWriter(fps=1000 // interval)
        anim.save(args.gif, writer=writer, dpi=160)
        print(f"Saved GIF: {args.gif}  (frames={len(frames)}, fps≈{1000//interval})")

    # optional MP4 (requires ffmpeg in PATH)
    if args.mp4:
        try:
            anim.save(args.mp4, writer="ffmpeg", fps=1000 // interval, dpi=160)
            print(f"Saved MP4: {args.mp4}")
        except Exception as e:
            print(f"[warn] MP4 export failed (need ffmpeg?): {e}")

if __name__ == "__main__":
    main()
