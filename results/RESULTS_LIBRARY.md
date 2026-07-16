# Results library (pointer)

Analysis **outputs** under `results/` are inventoried **outside** this FLITS tree:

```text
~/Data/Faber2026/results-library/INDEX.md
~/Data/Faber2026/results-library/_inventory/inventory.yaml
```

Override root with `FABER2026_RESULTS_LIBRARY`. Phase B may relocate bulky gitignored trees into the library and leave README stubs here; until then this file is a pointer only.

Catalog + builder live in the sibling Faber2026 manuscript repo:

```text
scripts/results_library_catalog.yaml
scripts/build_results_library_inventory.py
scripts/results_library.py   # results_slot(...)
```

Refresh from the Faber2026 checkout:

```bash
python3 scripts/build_results_library_inventory.py --dry-run
python3 scripts/build_results_library_inventory.py --link --force
```

See also [`../DATA_LOCATIONS.md`](../DATA_LOCATIONS.md) (Results library section) and [`../analysis/RESULTS_LIBRARY.md`](../analysis/RESULTS_LIBRARY.md).
