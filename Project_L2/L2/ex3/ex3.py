import os
from typing import Dict, List, Tuple

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


# ---------- I/O helpers ----------

def load_fasta_file(path: str) -> str:
    """
    Read a FASTA file (single or multi-record) and return the concatenated
    sequence as uppercase with whitespace removed. Header lines (>) are skipped.
    """
    chunks: List[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith(">"):
                continue
            chunks.append(line)
    return "".join(chunks).replace(" ", "").upper()


# ---------- analysis ----------

def window_freqs(
    seq: str,
    win: int = 30,
    alphabet: Tuple[str, ...] = ("A", "C", "G", "T"),
    denom_only_alphabet: bool = True,
) -> Tuple[List[int], Dict[str, List[float]]]:
    """
    Slide a window across `seq` and compute relative frequency for each symbol
    in `alphabet`. Returns 1-based window starts + per-symbol frequency lists.
    """
    n = len(seq)
    if win <= 0:
        raise ValueError("Window must be > 0.")
    if n < win:
        return [], {s: [] for s in alphabet}

    starts: List[int] = []
    freqs: Dict[str, List[float]] = {s: [] for s in alphabet}

    for i in range(0, n - win + 1):
        frag = seq[i:i + win]
        if denom_only_alphabet:
            denom = sum(1 for ch in frag if ch in alphabet)
        else:
            denom = win

        counts = {s: 0 for s in alphabet}
        for ch in frag:
            if ch in counts:
                counts[ch] += 1

        if denom == 0:
            rel = {s: 0.0 for s in alphabet}
        else:
            rel = {s: counts[s] / denom for s in alphabet}

        starts.append(i + 1)  # 1-based for user display
        for s in alphabet:
            freqs[s].append(rel[s])

    return starts, freqs


# ---------- GUI ----------

class NucFreqGUI(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("Sliding-Window Nucleotide Frequencies")
        self.geometry("1000x700")

        self._seq: str = ""
        self._path: str | None = None

        self._make_widgets()

    # UI layout
    def _make_widgets(self) -> None:
        bar = ttk.Frame(self)
        bar.pack(side=tk.TOP, fill=tk.X, padx=10, pady=8)

        ttk.Button(bar, text="Open FASTA…", command=self._on_open).grid(row=0, column=0, padx=5, sticky="w")

        self._file_var = tk.StringVar(value="No file selected")
        ttk.Label(bar, textvariable=self._file_var, width=60).grid(row=0, column=1, padx=5, sticky="w")

        ttk.Label(bar, text="Window:").grid(row=0, column=2, padx=(18, 5), sticky="e")
        self._win_var = tk.IntVar(value=30)
        ttk.Spinbox(bar, from_=1, to=10000, textvariable=self._win_var, width=8).grid(row=0, column=3, padx=5, sticky="w")

        self._denom_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            bar, text="Exclude non-ACGT from denominator", variable=self._denom_var
        ).grid(row=0, column=4, padx=(18, 5), sticky="w")

        ttk.Button(bar, text="Analyze & Plot", command=self._on_analyze).grid(row=0, column=5, padx=12, sticky="e")

        self._info = tk.StringVar(value="Load a FASTA file to begin.")
        ttk.Label(self, textvariable=self._info, foreground="#333").pack(side=tk.TOP, anchor="w", padx=12, pady=(0,8))

        # figure
        self._fig = Figure(figsize=(10, 5), dpi=100)
        self._ax = self._fig.add_subplot(111)
        self._prime_axes()

        self._canvas = FigureCanvasTkAgg(self._fig, master=self)
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        self._status = tk.StringVar(value="")
        ttk.Label(self, textvariable=self._status, anchor="w", relief=tk.SUNKEN).pack(side=tk.BOTTOM, fill=tk.X)

        # optional theme
        try:
            self.tk.call("source", "sun-valley.tcl")
            self.tk.call("set_theme", "light")
        except Exception:
            pass

    def _prime_axes(self) -> None:
        self._ax.clear()
        self._ax.set_title("Sliding-Window Relative Frequencies (A, C, G, T)")
        self._ax.set_xlabel("Window start position (1-based)")
        self._ax.set_ylabel("Relative frequency")
        self._ax.set_ylim(0, 1)
        self._ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.5)

    # actions
    def _on_open(self) -> None:
        p = filedialog.askopenfilename(
            title="Select FASTA file",
            filetypes=[("FASTA", "*.fa *.fasta *.fna *.ffn *.faa *.frn"), ("All files", "*.*")]
        )
        if not p:
            return
        try:
            seq = load_fasta_file(p)
            if not seq:
                raise ValueError("No bases found in file.")
            self._seq = seq
            self._path = p
            self._file_var.set(os.path.basename(p))
            alpha = "".join(sorted(set(self._seq)))
            self._info.set(f"Loaded sequence length: {len(self._seq)} | Alphabet: {alpha}")
            self._status.set("File loaded.")
        except Exception as exc:
            messagebox.showerror("FASTA error", str(exc))
            self._status.set("Failed to load file.")

    def _on_analyze(self) -> None:
        if not self._seq:
            messagebox.showwarning("No sequence", "Open a FASTA file first.")
            return

        try:
            w = int(self._win_var.get())
            if w <= 0:
                raise ValueError
        except Exception:
            messagebox.showerror("Invalid window", "Window must be a positive integer.")
            return

        alphabet = ("A", "C", "G", "T")
        pos, fr = window_freqs(
            self._seq,
            win=w,
            alphabet=alphabet,
            denom_only_alphabet=bool(self._denom_var.get()),
        )

        if not pos:
            self._info.set(f"Sequence length ({len(self._seq)}) is smaller than window size ({w}).")
            self._prime_axes()
            self._canvas.draw_idle()
            return

        self._info.set(
            f"Analyzed {len(pos)} windows | Window={w} | "
            f"Denominator={'ACGT only' if self._denom_var.get() else 'full window'}"
        )

        # draw
        self._prime_axes()
        for s in alphabet:
            self._ax.plot(pos, fr[s], label=s, linewidth=1.6)
        self._ax.legend(title="Nucleotide", ncol=len(alphabet), frameon=True)
        self._canvas.draw_idle()


def main() -> None:
    app = NucFreqGUI()
    app.mainloop()


if __name__ == "__main__":
    main()
