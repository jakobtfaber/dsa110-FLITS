# Machine inventory (CHIME–DSA codetections)

Canonical file: [`machine_inventory.yaml`](../../machine_inventory.yaml)

Initial inventory live-probed 2026-06-25 across jakob-mbp, iacobus, h17, h23,
hpcc, dsacamera, and arc (VOSpace), with scoped custody updates through
2026-07-21. Current access roles are defined in
[`DATA_LOCATIONS.md`](../../DATA_LOCATIONS.md).

**4-host migration plan:** [`MIGRATION_PLAN_4HOST.md`](MIGRATION_PLAN_4HOST.md) — Phase 2 design: [`PHASE2_DESIGN.md`](PHASE2_DESIGN.md) — log: [`MIGRATION_LOG.md`](MIGRATION_LOG.md)

## Quick reference

| Machine | migration_status | Role | Key path |
|---------|------------------|------|----------|
| **iacobus** | retired | drained staging source | `~/Research/_quarantine/CHIME_DSA_Codetections-drained-20260713/` |
| **jakob-mbp** | target | development + local replicas | `~/Developer/repos/github.com/jakobtfaber/dsa110-FLITS` |
| **arc** | target | fit-input store + CANFAR compute | `arc:home/jfaber/baseband_morphologies/chime_dsa_codetections` |
| **h17** | target | raw-data access + compute | `/data/Faber2026/data` |
| **h23** | retired | cold upstream (drain) | `/media/ubuntu/ssd/jfaber/chime_dsa_codetections` |
| **hpcc** | retired | Slurm batch (drain) | `/home/jfaber/flits/dsa110-FLITS` |
| **dsacamera** | retired | negligible codetection | — |

## Query with Python

```bash
# all paths on iacobus
python scripts/query_machine_inventory.py --machine iacobus

# every git repo entry
python scripts/query_machine_inventory.py --kind git_repo --json

# find DSA burst paths anywhere
python scripts/query_machine_inventory.py --path-contains DSA_bursts

# pending migration_map entries
python scripts/query_machine_inventory.py --migration-map --migration-status pending

# locations on retired hosts
python scripts/query_machine_inventory.py --migration-status retired

# everything targeting iacobus (locations + migration_map union)
python scripts/query_machine_inventory.py --migration-target iacobus

# retired-host coverage gate (exit 0 = all h23 codetection subtrees mapped)
python scripts/query_machine_inventory.py --check-retired-coverage
```

## Query with yq

```bash
yq '.migration' machine_inventory.yaml
yq '.migration_map[] | select(.status=="pending")' machine_inventory.yaml
yq '.canonical' machine_inventory.yaml
yq '.machines.iacobus.locations[] | select(.id=="research_chime_codetections")' machine_inventory.yaml
yq '.. | select(has("path")) | .path' machine_inventory.yaml | rg burst
```

## Refresh

Re-run SSH/`vls` probes and update `generated_utc` + affected subtrees. Companion manifest for iacobus sentinels: [`codetections_manifest.yaml`](../../codetections_manifest.yaml).
