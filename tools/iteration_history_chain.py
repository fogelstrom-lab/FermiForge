"""Join completed restart segments without changing either source history."""
from pathlib import Path
import math

META = ("local_iteration", "segment", "history_restart")


def read_history(path):
    lines = Path(path).read_text().splitlines()
    if not lines or not lines[0].startswith("#"):
        raise ValueError(f"missing history header: {path}")
    names = lines[0][1:].split()
    if len(set(names)) != len(names) or "iteration" not in names:
        raise ValueError(f"invalid history columns: {path}")
    rows = []
    for line in lines[1:]:
        if not line.strip() or line.startswith("#"):
            continue
        values = line.split()
        if len(values) != len(names) or not all(math.isfinite(float(v)) for v in values):
            raise ValueError(f"invalid history row: {path}")
        rows.append(dict(zip(names, values)))
    if not rows or [float(r["iteration"]) for r in rows] != list(range(1, len(rows)+1)):
        raise ValueError(f"history must contain consecutive iterations starting at 1: {path}")
    return names, rows


def concatenate_history(parent, current, output):
    pn, pr = read_history(parent)
    cn, cr = read_history(current)
    columns = [n for n in cn if n not in META]
    if columns != [n for n in pn if n not in META]:
        raise ValueError("parent and continuation history schemas differ")
    offset = len(pr)
    last_segment = int(float(pr[-1].get("segment", "1")))
    lines = ["# " + " ".join(columns + list(META))]
    for row in pr:
        lines.append(" ".join([row[n] for n in columns] + [
            row.get("local_iteration", row["iteration"]),
            row.get("segment", "1"), row.get("history_restart", "0")]))
    for i, row in enumerate(cr, 1):
        joined = dict(row, iteration=str(offset+i))
        lines.append(" ".join([joined[n] for n in columns] +
                              [str(i), str(last_segment+1), str(int(i==1))]))
    Path(output).write_text("\n".join(lines)+"\n")
    return offset, len(cr)
