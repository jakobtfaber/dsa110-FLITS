# Machine inventory (CHIME–DSA codetections)

Canonical file: [`machine_inventory.yaml`](../../machine_inventory.yaml)

The inventory records current project locations on jakob-mbp, h17, h23, hpcc,
dsacamera, and arc (VOSpace), with scoped custody updates through 2026-07-21.
Current access roles are defined in
[`DATA_LOCATIONS.md`](../../DATA_LOCATIONS.md).

## Quick reference

| Machine | migration_status | Role | Key path |
|---------|------------------|------|----------|
| **jakob-mbp** | target | development + local replicas | `~/Developer/repos/github.com/jakobtfaber/dsa110-FLITS` |
| **arc** | target | fit-input store + CANFAR compute | `arc:home/jfaber/baseband_morphologies/chime_dsa_codetections` |
| **h17** | target | raw-data access + compute | `/data/Faber2026/data` |
| **h23** | retired | cold upstream (drain) | `/media/ubuntu/ssd/jfaber/chime_dsa_codetections` |
| **hpcc** | retired | Slurm batch (drain) | `/home/jfaber/flits/dsa110-FLITS` |
| **dsacamera** | retired | negligible codetection | — |

## Query with Python

```bash
# every git repo entry
python scripts/query_machine_inventory.py --kind git_repo --json

# find DSA burst paths anywhere
python scripts/query_machine_inventory.py --path-contains DSA_bursts

# locations on retired hosts
python scripts/query_machine_inventory.py --migration-status retired
```

## Query with yq

```bash
yq '.canonical' machine_inventory.yaml
yq '.. | select(has("path")) | .path' machine_inventory.yaml | rg burst
```

## Refresh

Re-run SSH/`vls` probes and update `generated_utc` plus affected subtrees.
