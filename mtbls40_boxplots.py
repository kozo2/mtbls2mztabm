"""Plot per-metabolite abundance boxplots from the pickled mzTab-M model.

For each unique non-null InChI in small_molecule_summary, produce a boxplot
of the abundance_assay values grouped by metadata.study_variable.

Outputs: output/boxplots/<safe-inchi>.png
"""

from __future__ import annotations

import pickle
import re
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
PICKLE_PATH = ROOT / "output" / "MTBLS40.pkl"
OUT_DIR = ROOT / "output" / "boxplots"


def safe_name(inchi: str, max_len: int = 80) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", inchi)
    return s[:max_len].strip("_") or "metabolite"


def main() -> int:
    with PICKLE_PATH.open("rb") as fh:
        model = pickle.load(fh)

    study_vars = model.metadata.study_variable or []
    if not study_vars:
        print("No study_variable in metadata; cannot plot.")
        return 1

    # study_variable -> list of 0-based assay indices into abundance_assay
    sv_labels: List[str] = []
    sv_assay_idx: List[List[int]] = []
    for sv in study_vars:
        sv_labels.append(sv.name or f"sv[{sv.id}]")
        sv_assay_idx.append([a - 1 for a in (sv.assay_refs or [])])

    sml = model.small_molecule_summary or []
    by_inchi: Dict[str, List[Optional[float]]] = {}
    name_for: Dict[str, str] = {}
    for row in sml:
        inchi_list = row.inchi or []
        inchi = next((i for i in inchi_list if i and i.lower() != "null"), None)
        if not inchi:
            continue
        if inchi in by_inchi:
            # Skip duplicate rows for the same metabolite
            continue
        by_inchi[inchi] = list(row.abundance_assay or [])
        chem_name = (row.chemical_name or [""])[0] or "unknown"
        name_for[inchi] = chem_name

    if not by_inchi:
        print("No InChI-identified metabolites in small_molecule_summary.")
        return 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Plotting {len(by_inchi)} metabolites across {len(sv_labels)} study variables")

    for inchi, abundances in by_inchi.items():
        groups: List[List[float]] = []
        for assay_idx in sv_assay_idx:
            vals = [
                abundances[i]
                for i in assay_idx
                if 0 <= i < len(abundances) and abundances[i] is not None
            ]
            groups.append(vals)

        fig, ax = plt.subplots(figsize=(6, 4.5))
        ax.boxplot(groups, tick_labels=sv_labels, showfliers=True)
        ax.set_xlabel("study variable")
        ax.set_ylabel("abundance_assay")
        title = f"{name_for[inchi]}\n{inchi}"
        ax.set_title(title, fontsize=9)
        ax.grid(axis="y", linestyle=":", alpha=0.5)
        fig.tight_layout()

        out_path = OUT_DIR / f"{safe_name(name_for[inchi] + '_' + inchi)}.png"
        fig.savefig(out_path, dpi=120)
        plt.close(fig)

    print(f"Wrote {len(by_inchi)} boxplots to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
