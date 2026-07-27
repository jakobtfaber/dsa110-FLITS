#!/usr/bin/env python3
"""Query machine_inventory.yaml for paths, hosts, kinds, and migration state."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise SystemExit("PyYAML required: conda run -n py312 python -m pip install pyyaml") from exc

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INVENTORY = REPO_ROOT / "machine_inventory.yaml"

def load_inventory(path: Path) -> dict:
    with path.open(encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def walk_locations(node: dict, machine: str, hits: list[dict]) -> None:
    loc = dict(node)
    loc.setdefault("machine", machine)
    hits.append(loc)
    for child in loc.get("children") or []:
        if isinstance(child, str):
            hits.append({"machine": machine, "path": child})
        elif isinstance(child, dict):
            walk_locations(child, machine, hits)


def flatten(inventory: dict) -> list[dict]:
    rows: list[dict] = []
    for machine, meta in (inventory.get("machines") or {}).items():
        for loc in meta.get("locations") or []:
            walk_locations(loc, machine, rows)
    return rows


def filter_location_rows(inv: dict, args: argparse.Namespace) -> list[dict]:
    rows = flatten(inv)
    if args.machine:
        rows = [r for r in rows if r.get("machine") == args.machine]
    if args.kind:
        rows = [r for r in rows if r.get("kind") == args.kind]
    if args.path_contains:
        needle = args.path_contains.lower()
        rows = [r for r in rows if needle in str(r.get("path", "")).lower()]
    if args.migration_status:
        machines = inv.get("machines") or {}
        allowed = {
            name
            for name, meta in machines.items()
            if meta.get("migration_status") == args.migration_status
        }
        rows = [r for r in rows if r.get("machine") in allowed]
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--machine", help="Filter by machine id (h17, arc, ...)")
    parser.add_argument("--kind", help="Filter by location kind")
    parser.add_argument("--path-contains", help="Substring match on path")
    parser.add_argument(
        "--migration-status",
        help="Filter machines by migration_status",
    )
    parser.add_argument("--json", action="store_true", help="Emit JSON")
    args = parser.parse_args()

    inv = load_inventory(args.inventory)

    rows = filter_location_rows(inv, args)

    if args.json:
        json.dump(rows, sys.stdout, indent=2)
        sys.stdout.write("\n")
        return 0

    for row in rows:
        path = row.get("path", "?")
        kind = row.get("kind", "")
        size = row.get("size_human", "")
        extra = f" [{kind}]" if kind else ""
        size_s = f" {size}" if size else ""
        mig = row.get("migration_action")
        mig_s = f" ({mig}->{row.get('migration_target')})" if mig else ""
        prefix = row.get("machine", "")
        if str(path).startswith("arc:"):
            label = str(path)
        else:
            label = f"{prefix}:{path}"
        print(f"{label}{extra}{size_s}{mig_s}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
