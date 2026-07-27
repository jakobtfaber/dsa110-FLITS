# 4-host migration plan (CHIME–DSA codetections)

**Status:** Phase 5 complete 2026-06-25 (retired hosts quarantined) — log: [`MIGRATION_LOG.md`](MIGRATION_LOG.md)  
**Inventory:** [`machine_inventory.yaml`](../../machine_inventory.yaml) (`migration_map`, per-path `migration_*` fields)  
**Sentinels:** [`codetections_manifest.yaml`](../../codetections_manifest.yaml)

Consolidate code, data, docs, and run artifacts onto **four hosts only**. Retire h23, hpcc, and dsacamera as *homes* (move-only quarantine — never delete).

## Target hosts

| Host | Role | Primary paths |
|------|------|---------------|
| **jakob-mbp** | Code + manuscripts + small local replicas | `~/Developer/repos/github.com/jakobtfaber/dsa110-FLITS` |
| **iacobus** | Then-current staging tree + iCloud uploader | `~/Research/CHIME_DSA_Codetections` (218G) |
| **arc** | Institutional burst `.npy`, baseband, CANFAR compute | `arc:home/jfaber/baseband_morphologies/chime_dsa_codetections` |
| **h17** | OVRO compute, docker CHIME, staging archives | `/data/research/.../chime-dsa-codetections` (canonical; includes arc archive + upchan) |

**Git source of truth:** GitHub `jakobtfaber/dsa110-FLITS` — not hpcc, arc checkout, or h23 trees.

## Retired hosts (drain → quarantine)

| Host | Action after drain |
|------|-------------------|
| **h23** | Read-only; unique bytes → iacobus/arc; then `mv` to `_quarantine/h23-drain-*` |
| **hpcc** | Pull run artifacts → git/iacobus; quarantine `/home/jfaber/flits/` |
| **dsacamera** | No codetection content; mark decommissioned in inventory |

**Out of scope:** nihari, Faber2024_bursts, h23 `/dataz/dsa110/T3` (59T raw pipeline).

## Data tiering

| Asset class | Primary | Mirror / compute |
|-------------|---------|------------------|
| FLITS code | jakob-mbp + GitHub | h17 optional clone for docker workflows |
| Burst pickles (61G) | iacobus | iCloud placeholders on jakob-mbp |
| DSA/CHIME `.npy` for fits | arc `DSA_bursts` / `CHIME_bursts` | jakob-mbp local replica ≤5G; iacobus optional |
| Waterfalls, npys, scattering | iacobus | merge hpcc JSON into scattering_results |
| Archive (131G+) | iacobus | dedupe vs arc `OLD_CHIME_DSA_Codetections` (77G) |
| Baseband / singlebeam | h17 compute workspace | promote finished products → arc |
| arc trash archive | h17 `arc_archive_2026-06` | optional copy → iacobus after dedupe |
| CANFAR `flits_run` | arc ephemeral | regenerate; do not migrate |

**jakob-mbp constraint:** 25G free (98% disk) — never materialize the 218G iacobus tree locally.

## Phases

### Phase 0 — Manifest (1 day)

1. Refresh `machine_inventory.yaml` and `codetections_manifest.yaml`.
2. Every path tagged in `migration_map` with `action`, `target_host`, `phase`, `status`.
3. Disk budget: jakob-mbp, arc quota (~200G), iacobus free space.

**Exit:** no unmapped codetection paths on retired hosts.

### Phase 1 — Code → jakob-mbp + GitHub

| Source | Action |
|--------|--------|
| jakob-mbp FLITS @ `7e4c0c97` | Keep canonical |
| hpcc `/home/jfaber/flits/dsa110-FLITS` | Pull `_a1_fits/*.json` (35) → repo or iacobus; quarantine hpcc tree |
| arc `arc:home/jfaber/dsa110-FLITS` | Diff vs GitHub; rescue CANFAR-only edits via PR; stop dev checkout |
| arc `chime_dsa_codetections/FLITS` (~5G) | Diff vs GitHub; dedupe; not canonical |
| h23 `chime_dsa_codetections/` scripts | Worth keeping → FLITS; rest quarantine |
| Legacy Documents clones | Drain → `Developer/repos/` + `_quarantine/` |

### Phase 2 — Data drain (iacobus as uploader)

**Design doc:** [`PHASE2_DESIGN.md`](PHASE2_DESIGN.md) — audit-first; legacy paths corrected; wave order; burst_npys + OLD_CHIME policies.

Use iacobus (not jakob-mbp) for h23→iacobus rsync (`scripts/migration/h23_to_iacobus.sh`, adapted from legacy `h23_resilient_transfer.sh`).

| h23 source | Target |
|------------|--------|
| `.../chime_dsa_codetections/bursts` (7.5G) | iacobus `burst_npys` (reconcile 32 vs 52 files first) |
| `.../data/stokes_I_npys` (983M) | iacobus `dsa_fullstokes_waterfalls` |
| `.../dm_budget` (329M) | iacobus `dm_budget` |
| `.../scattering` (24M) | iacobus `scattering_results` |
| `.../dm` (732K) | iacobus `dm_budget` (merge) |
| `.../localizations` (8K) | iacobus `metadata` |
| `.../archive/OLD_CHIME_DSA_Codetections` (79G) | iacobus `archive/` (verify sentinels) |
| `.../burstprop_paper/bursts` (42G) | iacobus `archive/burstprop_paper` |
| `.../archive/dsa110-scat` (5.9G) | git repo + iacobus if needed |

### Phase 3 — arc ↔ iacobus dedupe

**Design doc:** [`PHASE3_DESIGN.md`](PHASE3_DESIGN.md) — audit-first; quota guard; OLD_CHIME sentinel policy; no bulk 218G upload.

**Started:** 2026-06-25 (parallel with Phase 2 rsync). Initial audit: [`reports/phase3_audit.json`](../../reports/phase3_audit.json).

| Audit id | arc | iacobus / jakob | Initial recommendation |
|----------|-----|-----------------|------------------------|
| `arc_dsa_bursts` | 2.9G | jakob local 2.8G | sync local replica gaps |
| `arc_chime_bursts` | 5.9G | iacobus burst_npys 21G | reconcile namespaces |
| `arc_old_chime_dedupe` | 77G | iacobus 105G | **iacobus canonical** (sentinel pass) |
| `arc_flits_checkout` | 211M | jakob/GitHub | diff vs GitHub |
| `arc_codetection_flits_tree` | 5G | jakob/GitHub | diff vs GitHub |

1. Path-level diff: `DSA_bursts`, `CHIME_bursts` vs iacobus/jakob local (`audit_arc_delta.py`).
2. Pick canonical for `OLD_CHIME_DSA_Codetections` (arc 77G vs iacobus 105G) by sentinel hash + layout map.
3. Do **not** bulk-upload 218G to arc without quota check (~200G cap).
4. Phase 1 dedupe ids `arc_flits_checkout`, `arc_codetection_flits_tree` — audit in Phase 3.1.

### Phase 4 — h17 (already a target)

**Design doc:** [`PHASE4_DESIGN.md`](PHASE4_DESIGN.md) — audit-first; keep compute; optional arc-trash copy; no h23/rsync collision.

**Started:** 2026-06-25 (parallel with Phase 2 rsync / Phase 3 dedupe). Initial audit: [`reports/phase4_audit.json`](../../reports/phase4_audit.json).

| Audit id | h17 | iacobus compare | Initial recommendation |
|----------|-----|-----------------|------------------------|
| `h17_compute_workspace` | 29G / 339 f | — | **keep** |
| `h17_upchan_products` | 473M / 11 f | — | **keep**; promote when stable |
| `h17_arc_archive_copy` | 36G / 1924 f | target missing | **dedupe_then_copy** after Phase 3 |
| `h17_ubuntu_stub` | 0 f | — | **remove_stub** |
| `iacobus_chime_canfar_archive` | — | 2.7G vs research archive | **dedupe_into_research** |

1. Keep compute workspace (29G) and upchan products (473M).
2. Promote finished artifacts → arc or iacobus only when stable.
3. Optional: copy `arc_archive_2026-06` → iacobus after OLD_CHIME dedupe (`h17_arc_archive_copy`).
4. Remove empty stubs (`/data/ubuntu/chime-dsa-codetections`, `/data/jfaber/chime_singlebeam`).
5. **Do not** invoke Phase 2 h23→iacobus transfers from Phase 4 work.

### Phase 5 — Quarantine + decommission

1. `_quarantine/README.md` entry per retired path with one-line restore.
2. Update `DATA_LOCATIONS.md` to 4-host model only.
3. Inventory query: zero `status: pending` on retired-host `migration_map` entries.

## Verification gates

Per subtree:

- File count + size match manifest
- SHA-256 sentinel from `codetections_manifest.yaml`
- FLITS configs point at arc paths or documented iacobus mirror
- `python scripts/query_machine_inventory.py --migration-status pending` returns empty
- `python scripts/query_machine_inventory.py --check-retired-coverage` exits 0 (all retired-host codetection paths covered by `migration_map`)

## Query migration state

```bash
# all pending migrations
python scripts/query_machine_inventory.py --migration-status pending

# everything targeting iacobus (locations + migration_map)
python scripts/query_machine_inventory.py --migration-target iacobus
python scripts/query_machine_inventory.py --migration-map --migration-target iacobus

# retired-host coverage gate
python scripts/query_machine_inventory.py --check-retired-coverage

yq '.migration_map[] | select(.status=="pending")' machine_inventory.yaml
yq '.migration.target_hosts' machine_inventory.yaml
```

## Ongoing ops (post-migration)

| Need | Host |
|------|------|
| Edit FLITS | jakob-mbp → push GitHub |
| Bulk storage / iCloud | iacobus |
| Pipeline `.npy` / CANFAR GPU fits | arc |
| CHIME docker / upchan / 2×2080 Ti | h17 |
| Ad-hoc Slurm | hpcc burst only — no persistent project data |

## Risks

- jakob-mbp disk full — block local materialization of iCloud tree
- arc 200G quota — dedupe before upload
- burst_npys count mismatch (32 vs 52) — reconcile before h23 quarantine
- hpcc stale git SHA — artifact source only, not a branch
