"""Read output/MTBLS40.mztab via pymzTab-m and pickle the model.

Also reloads the pickle to confirm it can be deserialized.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import mztab_m_io as mztabm
from mztab_m_io.model.validation import MessageType

ROOT = Path(__file__).resolve().parent
MZTAB_PATH = ROOT / "output" / "MTBLS40.mztab"
PICKLE_PATH = ROOT / "output" / "MTBLS40.pkl"


def main() -> int:
    result = mztabm.read(str(MZTAB_PATH), format="tsv", auto_complete_ids=True)
    errors = [m for m in result.messages if m.message_type == MessageType.ERROR]
    print(f"Read {MZTAB_PATH}: success={result.success}, errors={len(errors)}")
    if not result.mztabm:
        for m in errors[:10]:
            print(f"  [ERROR] {m.message}")
        return 1

    PICKLE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with PICKLE_PATH.open("wb") as fh:
        pickle.dump(result.mztabm, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Wrote pickle: {PICKLE_PATH} ({PICKLE_PATH.stat().st_size:,} bytes)")

    with PICKLE_PATH.open("rb") as fh:
        loaded = pickle.load(fh)

    print(f"Reloaded type: {type(loaded).__module__}.{type(loaded).__name__}")
    print(f"  mzTab-ID: {loaded.metadata.mztab_id}")
    print(f"  title:    {loaded.metadata.title}")
    print(f"  samples:  {len(loaded.metadata.sample or [])}")
    print(f"  ms_runs:  {len(loaded.metadata.ms_run or [])}")
    print(f"  assays:   {len(loaded.metadata.assay or [])}")
    print(f"  study_variables: {len(loaded.metadata.study_variable or [])}")
    print(f"  SML rows: {len(loaded.small_molecule_summary or [])}")
    print(f"  SMF rows: {len(loaded.small_molecule_feature or [])}")
    print(f"  SME rows: {len(loaded.small_molecule_evidence or [])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
