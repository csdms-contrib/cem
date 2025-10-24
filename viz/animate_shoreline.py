#!/usr/bin/env python3
import argparse, glob, io, re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from pathlib import Path

def load_shoreline(path):
    txt = Path(path).read_text(encoding="utf-8", errors="ignore").replace("%","").strip()
    txt = re.sub(r"\s+", " ", txt)
    arr = np.loadtxt(io.StringIO(txt))
    return arr if arr.ndim == 1 else arr[:,1] if arr.shape[1] > 1 else arr[:,0]

def main():
    p = argparse.ArgumentParser(description="Animate CEM shoreline snapshots.")
    p.add_argument("--input", default="ShorePos_*.dat", help="Glob for shoreline files")
    p.add_argument("--ds", type=float, default=None, help="Alongshore spacing (m). If not set, x is cell index.")
    p.add_argument("--gif", default="shoreline.gif", help="Output GIF filename")
    p.add_argument("--interval", type=int, default=300, help="Frame interval (ms)")
    args = p.parse_args()

    files = sorted(glob.glob(args.input))
    if not files:
        raise SystemExit(f"No files match: {args.input}")

    ys = [load_shoreline(f) for f in files]
    n = len(ys[0])
    x = np.arange(n) if args.ds is None else np.arange(n) * args.ds

    ymin = min(map(np.min, ys))
    ymax = max(map(np.max, ys))

    fig, ax = plt.subplots()
    (line,) = ax.plot(x, ys[0])
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("Alongshore distance (m)" if args.ds is not None else "Alongshore cell index")
    ax.set_ylabel("Shoreline position y")
    title = ax.set_title(Path(files[0]).stem)

    def update(i):
        line.set_ydata(ys[i])
        title.set_text(Path(files[i]).stem)
        return (line, title)

    anim = FuncAnimation(fig, update, frames=len(ys), interval=args.interval, blit=True)
    anim.save(args.gif, writer=PillowWriter(fps=1000//args.interval))
    print(f"Saved {args.gif}")

if __name__ == "__main__":
    main()
