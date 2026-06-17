# CEM terminal session — 2025-10-24

- host: waveletswave-3
- shell: zsh
- cwd: /Users/benthosyy/Desktop/1024CEM/cem/_build/runs/more

## Commands

```bash
 1001  24.10.2025 17:17  git clone https://github.com/waveletswave/cem.git
 1002  24.10.2025 17:18  code cem
 1003  24.10.2025 17:25  mkdir _build
 1004  24.10.2025 17:26  cd _build
 1005  24.10.2025 17:26  cmake ..
 1006  24.10.2025 17:26  make -j
 1007  24.10.2025 17:27  cd ~/Desktop/1024CEM
 1008  24.10.2025 17:27  ls -la
 1009  24.10.2025 17:28  cd ~/Desktop/1024CEM/cem
 1010  24.10.2025 17:28  ls -la | grep -E 'CMakeLists.txt|main.c'
 1011  24.10.2025 17:29  cd ~/Desktop/1024CEM/cem
 1012  24.10.2025 17:29  cmake -S . -B _build\n
 1013  24.10.2025 17:31  cd ~/Desktop/1024CEM/cem
 1014  24.10.2025 17:31  cmake --build _build -j\n
 1015  24.10.2025 17:39  cd ~/Desktop/1024CEM/cem
 1016  24.10.2025 17:39  rm -rf _build\n
 1017  24.10.2025 17:39  cmake -S . -B _build
 1018  24.10.2025 17:40  cmake --build _build -j
 1019  24.10.2025 17:41  cd ~/Desktop/1024CEM/cem
 1020  24.10.2025 17:41  cmake --build _build -- -j1
 1021  24.10.2025 17:43  cd ~/Desktop/1024CEM/cem/_build
 1022  24.10.2025 17:43  mkdir -p CMakeFiles/bmi_cem.dir/cem CMakeFiles/bmi_cem.dir/bmi \\n         CMakeFiles/bmi_cem-static.dir/cem CMakeFiles/bmi_cem-static.dir/bmi
 1023  24.10.2025 17:43  cmake --build . -j
 1024  24.10.2025 17:47  cd ~/Desktop/1024CEM/cem/_build
 1025  24.10.2025 17:47  mkdir -p runs/minimal
 1026  24.10.2025 17:47  cd runs/minimal
 1027  24.10.2025 17:47  ../../cem
 1028  24.10.2025 17:49  ls -lh
 1029  24.10.2025 17:50  head -n 5 ShorePos_1.dat
 1030  24.10.2025 17:50  head -n 5 CEM_365.out
 1031  24.10.2025 17:51  cd ~/Desktop/1024CEM/cem/_build/runs/minimal\n\npython3 - <<'PY'\nimport glob, io, re, numpy as np, matplotlib.pyplot as plt\ndef load_shoreline(path):\n    # Read as text, strip any '%' and compress whitespace, then load as float array\n    txt = open(path, 'r', encoding='utf-8', errors='ignore').read()\n    txt = txt.replace('%','').strip()\n    buf = io.StringIO(re.sub(r'\s+', ' ', txt))\n    arr = np.loadtxt(buf)         # 1D array of shoreline positions y(i)\n    return arr\n\nfiles = sorted(glob.glob('ShorePos_*.dat'))\nif not files:\n    raise SystemExit("No ShorePos_*.dat files found.")\n\nplt.figure()\nfor f in files:\n    y = load_shoreline(f)\n    x = np.arange(len(y))         # alongshore cell index; replace with meters if you know ds\n    plt.plot(x, y, label=f)\n\nplt.xlabel("Alongshore cell index")\nplt.ylabel("Shoreline position y")\nplt.title("CEM Shoreline Snapshots")\nplt.legend()\nplt.tight_layout()\nplt.savefig("shoreline_snapshots.png", dpi=220)\nprint("Saved shoreline_snapshots.png")\nPY\n
 1032  24.10.2025 17:56  cd ~/Desktop/1024CEM/cem
 1033  24.10.2025 17:56  mkdir -p viz
 1034  24.10.2025 17:56  cat > viz/requirements.txt <<'TXT'\n\ncat > viz/requirements.txt <<'TXT'\nnumpy\nmatplotlib\npillow   # for GIF export\nTXT\n\n# (optional) install into your environment\npython3 -m pip install -r viz/requirements.txt
 1035  24.10.2025 17:58  cd ~/Desktop/1024CEM/cem\nrm -f viz/requirements.txt\ncat <<'TXT' > viz/requirements.txt\nnumpy\nmatplotlib\npillow   # for GIF export\nTXT\n
 1036  24.10.2025 17:58  python3 -m pip install -r viz/requirements.txt
 1037  24.10.2025 17:59  cat <<'PY' > viz/plot_shoreline.py\n#!/usr/bin/env python3\nimport argparse, glob, io, re\nimport numpy as np\nimport matplotlib.pyplot as plt\nfrom pathlib import Path\n\ndef load_shoreline(path, cols):\n    txt = Path(path).read_text(encoding="utf-8", errors="ignore")\n    txt = txt.replace("%","").strip()\n    txt = re.sub(r"\s+", " ", txt)\n    arr = np.loadtxt(io.StringIO(txt))\n    if arr.ndim == 1:\n        y = arr\n    else:\n        y = arr[:, cols[1]] if arr.shape[1] > 1 else arr[:,0]\n    return y\n\ndef main():\n    p = argparse.ArgumentParser(description="Plot CEM shoreline snapshots.")\n    p.add_argument("--input", default="ShorePos_*.dat", help="Glob for shoreline files")\n    p.add_argument("--cols", type=int, nargs=2, default=[0,1], help="Zero-based columns for x,y (used if file has multiple columns)")\n    p.add_argument("--ds", type=float, default=None, help="Alongshore spacing (m). If not set, x is cell index.")\n    p.add_argument("--labels", default=None, help="Comma-separated legend labels (same order as sorted files)")\n    p.add_argument("--out", default="shoreline_snapshots.png", help="Output PNG filename")\n    args = p.parse_args()\n\n    files = sorted(glob.glob(args.input))\n    if not files:\n        raise SystemExit(f"No files match: {args.input}")\n\n    labels = args.labels.split(",") if args.labels else None\n    if labels and len(labels) != len(files):\n        raise SystemExit("Number of labels must match number of files.")\n\n    xs = None\n    plt.figure()\n    for i, f in enumerate(files):\n        y = load_shoreline(f, args.cols)\n        if xs is None:\n            n = len(y)\n            xs = np.arange(n) if args.ds is None else np.arange(n) * args.ds\n        lab = labels[i] if labels else Path(f).stem\n        plt.plot(xs, y, label=lab)\n\n    plt.xlabel("Alongshore distance (m)" if args.ds is not None else "Alongshore cell index")\n    plt.ylabel("Shoreline position y")\n    plt.title("CEM Shoreline Snapshots")\n    plt.legend()\n    plt.tight_layout()\n    plt.savefig(args.out, dpi=220)\n    print(f"Saved {args.out}")\n\nif __name__ == "__main__":\n    main()\nPY\nchmod +x viz/plot_shoreline.py\n
 1038  24.10.2025 17:59  cat <<'PY' > viz/animate_shoreline.py\n#!/usr/bin/env python3\nimport argparse, glob, io, re\nimport numpy as np\nimport matplotlib.pyplot as plt\nfrom matplotlib.animation import FuncAnimation, PillowWriter\nfrom pathlib import Path\n\ndef load_shoreline(path):\n    txt = Path(path).read_text(encoding="utf-8", errors="ignore").replace("%","").strip()\n    txt = re.sub(r"\s+", " ", txt)\n    arr = np.loadtxt(io.StringIO(txt))\n    return arr if arr.ndim == 1 else arr[:,1] if arr.shape[1] > 1 else arr[:,0]\n\ndef main():\n    p = argparse.ArgumentParser(description="Animate CEM shoreline snapshots.")\n    p.add_argument("--input", default="ShorePos_*.dat", help="Glob for shoreline files")\n    p.add_argument("--ds", type=float, default=None, help="Alongshore spacing (m). If not set, x is cell index.")\n    p.add_argument("--gif", default="shoreline.gif", help="Output GIF filename")\n    p.add_argument("--interval", type=int, default=300, help="Frame interval (ms)")\n    args = p.parse_args()\n\n    files = sorted(glob.glob(args.input))\n    if not files:\n        raise SystemExit(f"No files match: {args.input}")\n\n    ys = [load_shoreline(f) for f in files]\n    n = len(ys[0])\n    x = np.arange(n) if args.ds is None else np.arange(n) * args.ds\n\n    ymin = min(map(np.min, ys))\n    ymax = max(map(np.max, ys))\n\n    fig, ax = plt.subplots()\n    (line,) = ax.plot(x, ys[0])\n    ax.set_xlim(x.min(), x.max())\n    ax.set_ylim(ymin, ymax)\n    ax.set_xlabel("Alongshore distance (m)" if args.ds is not None else "Alongshore cell index")\n    ax.set_ylabel("Shoreline position y")\n    title = ax.set_title(Path(files[0]).stem)\n\n    def update(i):\n        line.set_ydata(ys[i])\n        title.set_text(Path(files[i]).stem)\n        return (line, title)\n\n    anim = FuncAnimation(fig, update, frames=len(ys), interval=args.interval, blit=True)\n    anim.save(args.gif, writer=PillowWriter(fps=1000//args.interval))\n    print(f"Saved {args.gif}")\n\nif __name__ == "__main__":\n    main()\nPY\nchmod +x viz/animate_shoreline.py\n
 1039  24.10.2025 18:01  cd ~/Desktop/1024CEM/cem/_build/runs/minimal
 1040  24.10.2025 18:01  python3 ../../viz/plot_shoreline.py --input "ShorePos_*.dat" --out shoreline_snapshots.png
 1041  24.10.2025 18:02  python3 ~/Desktop/1024CEM/cem/viz/plot_shoreline.py --input "ShorePos_*.dat" --out shoreline_snapshots.png
 1042  24.10.2025 18:03  python3 ~/Desktop/1024CEM/cem/viz/animate_shoreline.py --input "ShorePos_*.dat" --gif shoreline.gif
 1043  24.10.2025 18:10  cd ~/Desktop/1024CEM/cem
 1044  24.10.2025 18:10  grep -Rn "ShorePos_" -n .
 1045  24.10.2025 18:21  cd ~/Desktop/1024CEM/cem
 1046  24.10.2025 18:21  grep -Rn "SAVE_SPACING\|SaveSpacing" cem
 1047  24.10.2025 18:36  cd ~/Desktop/1024CEM/cem
 1048  24.10.2025 18:36  cmake --build _build -j
 1049  24.10.2025 18:36  mkdir -p _build/runs/more
 1050  24.10.2025 18:37  cd _build/runs/more
 1051  24.10.2025 18:37  ../../cem
 1052  24.10.2025 18:37  ls ShorePos_*.dat | wc -l
 1053  24.10.2025 18:39  python3 ~/Desktop/1024CEM/cem/viz/plot_shoreline.py \\n  --input "ShorePos_*.dat" \\n  --out shoreline_snapshots.png\nopen shoreline_snapshots.png
 1054  24.10.2025 18:39  python3 ~/Desktop/1024CEM/cem/viz/animate_shoreline.py \\n  --input "ShorePos_*.dat" \\n  --gif shoreline_true.gif \\n  --fps 6\nopen shoreline_true.gif
 1055  24.10.2025 18:42  python3 - <<'PY'\nimport glob, numpy as np, matplotlib.pyplot as plt, pathlib, re, io\ndef load(path):\n    txt = open(path,'r',encoding='utf-8',errors='ignore').read().replace('%','').strip()\n    arr = np.loadtxt(io.StringIO(re.sub(r'\s+',' ',txt)))\n    return arr if arr.ndim==1 else arr[:,1] if arr.shape[1]>1 else arr[:,0]\n\nfiles = sorted(glob.glob("ShorePos_*.dat"),\n               key=lambda s:int(re.search(r'(\d+)', pathlib.Path(s).stem).group(1)))\nk = 8                                  \nsel_idx = np.linspace(0, len(files)-1, k).astype(int)\nsel = [files[i] for i in sel_idx]\n\nxs = None\nplt.figure(figsize=(9,4))\nfor f in sel:\n    y = load(f)\n    if xs is None: xs = np.arange(len(y))\n    plt.plot(xs, y, alpha=0.7, label=pathlib.Path(f).stem)\nplt.xlabel("Alongshore cell index"); plt.ylabel("Shoreline position y")\nplt.title("CEM Shoreline Snapshots (subset)")\nplt.legend(ncol=2, fontsize=8, framealpha=0.8)\nplt.tight_layout(); plt.savefig("shoreline_subset.png", dpi=220)\nprint("Saved shoreline_subset.png with", len(sel), "curves")\nPY\nopen shoreline_subset.png
 1056  24.10.2025 18:43  python3 - <<'PY'\nimport glob, numpy as np, matplotlib.pyplot as plt, re, io, pathlib\ndef load(path):\n    txt=open(path,'r',encoding='utf-8',errors='ignore').read().replace('%','').strip()\n    import re, io; txt=re.sub(r'\s+',' ',txt)\n    a=np.loadtxt(io.StringIO(txt)); return a if a.ndim==1 else a[:,1] if a.shape[1]>1 else a[:,0]\n\nfiles=sorted(glob.glob("ShorePos_*.dat"),\n             key=lambda s:int(re.search(r'(\d+)', pathlib.Path(s).stem).group(1)))\nYs=[load(f) for f in files]\nYs=np.vstack(Ys)\nx=np.arange(Ys.shape[1])\nmu=Ys.mean(axis=0); sd=Ys.std(axis=0)\n\nplt.figure(figsize=(9,4))\nplt.fill_between(x, mu-sd, mu+sd, alpha=0.25, label='mean ± 1σ')\nplt.plot(x, mu, lw=2, label='mean')\nplt.plot(x, Ys[0],  lw=1.2, label=pathlib.Path(files[0]).stem)\nplt.plot(x, Ys[-1], lw=1.2, label=pathlib.Path(files[-1]).stem)\nplt.xlabel("Alongshore cell index"); plt.ylabel("y")\nplt.title("CEM shoreline: mean ± std + first/last")\nplt.legend(framealpha=0.8); plt.tight_layout(); plt.savefig("shoreline_band.png", dpi=220)\nprint("Saved shoreline_band.png")\nPY\nopen shoreline_band.png\n
 1057  24.10.2025 18:43  python3 - <<'PY'\nimport glob, numpy as np, matplotlib.pyplot as plt, re, io, pathlib\ndef load(path):\n    txt=open(path,'r',encoding='utf-8',errors='ignore').read().replace('%','').strip()\n    import re, io; txt=re.sub(r'\s+',' ',txt)\n    a=np.loadtxt(io.StringIO(txt)); return a if a.ndim==1 else a[:,1] if a.shape[1]>1 else a[:,0]\nfiles=sorted(glob.glob("ShorePos_*.dat"),\n             key=lambda s:int(re.search(r'(\d+)', pathlib.Path(s).stem).group(1)))\nY=np.vstack([load(f) for f in files])           # shape: [time, space]\nplt.figure(figsize=(8,4))\nim=plt.imshow(Y, aspect='auto', origin='lower', interpolation='nearest')\nplt.colorbar(im, label='Shoreline y')\nplt.ylabel('Snapshot index (time)'); plt.xlabel('Alongshore cell index')\nplt.title('CEM shoreline evolution (heatmap)')\nplt.tight_layout(); plt.savefig('shoreline_heatmap.png', dpi=220)\nprint("Saved shoreline_heatmap.png")\nPY\nopen shoreline_heatmap.png\n
 1058  24.10.2025 18:44  python3 - <<'PY'\nimport glob, numpy as np, matplotlib.pyplot as plt, re, io, pathlib\ndef load(p):\n    t=open(p,'r',encoding='utf-8',errors='ignore').read().replace('%','').strip()\n    import re, io; t=re.sub(r'\s+',' ',t)\n    a=np.loadtxt(io.StringIO(t)); return a if a.ndim==1 else a[:,1] if a.shape[1]>1 else a[:,0]\nfiles=sorted(glob.glob("ShorePos_*.dat"),\n             key=lambda s:int(re.search(r'(\d+)', pathlib.Path(s).stem).group(1)))\nY=[load(f) for f in files]\nx=np.arange(len(Y[0]))\nbase=Y[0]\nsel_idx=np.linspace(0,len(Y)-1,6).astype(int)   \nplt.figure(figsize=(9,4))\nfor i in sel_idx:\n    plt.plot(x, Y[i]-base, label=pathlib.Path(files[i]).stem, alpha=0.8)\nplt.axhline(0, ls='--', lw=1)\nplt.xlabel("Alongshore cell index"); plt.ylabel("Δy (vs first)")\nplt.title("CEM shoreline change relative to first snapshot")\nplt.legend(ncol=2, fontsize=8, framealpha=0.8)\nplt.tight_layout(); plt.savefig("shoreline_delta_subset.png", dpi=220)\nprint("Saved shoreline_delta_subset.png")\nPY\nopen shoreline_delta_subset.png
 1059  24.10.2025 18:48  python3 ~/Desktop/1024CEM/cem/viz/plot_shoreline.py \\n  --input "ShorePos_*.dat" --outdir . --ds 100 --subset 8 --legend\n
 1060  24.10.2025 18:51  python3 ~/Desktop/1024CEM/cem/viz/animate_shoreline_enhanced.py \\n  --input "ShorePos_*.dat" --gif shoreline_true.gif --interval 400
 1061  24.10.2025 18:51  python3 ~/Desktop/1024CEM/cem/viz/animate_shoreline_enhanced.py \\n  --input "ShorePos_*.dat" --gif s.gif --mp4 s.mp4 --fps 6 --ylim 0 20
 1062  24.10.2025 18:52  python3 ~/Desktop/1024CEM/cem/viz/plot_shoreline_enhanced.py \\n  --input "ShorePos_*.dat" --outdir . --ds 100 --subset 8 --legend\n
 1063  24.10.2025 18:54  python3 ~/Desktop/1024CEM/cem/viz/animate_shoreline_enhanced.py \\n  --input "ShorePos_*.dat" --gif shoreline_true.gif --interval 400
 1064  24.10.2025 19:07  fc -ln 1 > ~/Desktop/cem_cmds_$(date +%Y%m%d).txt
 1065  24.10.2025 19:09  echo 'setopt EXTENDED_HISTORY' >> ~/.zshrc\nsource ~/.zshrc
 1066  24.10.2025 19:10  HISTTIMEFORMAT='%Y-%m-%d %H:%M:%S  '  # pretty-print when using history\nbuiltin history -E 1 > ~/Desktop/cem_cmds_$(date +%Y%m%d)_ts.txt
 1067  24.10.2025 19:10  HISTTIMEFORMAT='%Y-%m-%d %H:%M:%S  ' history -E 1 \\n| awk -v d="$(date +%Y-%m-%d)" '$0 ~ d' \\n> ~/Desktop/cem_cmds_today_$(date +%Y%m%d).txt
 1068  24.10.2025 19:11  echo 'setopt EXTENDED_HISTORY' >> ~/.zshrc\nsource ~/.zshrc
 1069  24.10.2025 19:11  HISTTIMEFORMAT='%Y-%m-%d %H:%M:%S  ' history -E 1 \\n| awk -v d="$(date +%Y-%m-%d)" '$0 ~ d' \\n> ~/Desktop/cem_cmds_today_$(date +%Y%m%d).txt
```
