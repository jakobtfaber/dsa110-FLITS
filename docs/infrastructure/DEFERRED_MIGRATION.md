# Deferred migration work (post Phase 5)

**Status:** Phases 1–5 closed on `main` @ `a8c2b004` (PR #67).  
**Inventory gate:** `python scripts/query_machine_inventory.py --migration-status pending` → empty.  
**Policy:** move-only; no bulk transfers without explicit approval. See [`MIGRATION_PLAN_4HOST.md`](MIGRATION_PLAN_4HOST.md).

The items below record work that was initially deferred. D5 later completed;
remaining skipped items are not blockers for the current host model.

---

## D1 — `CHIME_bursts` cross-namespace reconcile

| Side | Path | Size (audit 2026-06-25) |
|------|------|-------------------------|
| arc (fits) | `arc:…/data/CHIME_bursts` | 60 f / 6.3 G |
| iacobus (archive) | `~/Research/CHIME_DSA_Codetections/burst_npys` | 218 f / 30.7 G |

**Finding:** arc holds fit-ready `.npy` under `dmphase/` and `dmtransform/` namespaces (24 codetection basenames in [`reports/phase3_chime_basename_inventory.csv`](../../reports/phase3_chime_basename_inventory.csv)); iacobus `burst_npys` uses a mixed nickname/TNS namespace — **zero basename overlap** with arc inventory rows.

**Map generated 2026-06-26** (read-only; no data movement):

```bash
python scripts/migration/map_chime_bursts_namespaces.py --stdout
```

Artifacts: [`reports/d1_chime_burst_map.csv`](../../reports/d1_chime_burst_map.csv), [`reports/d1_chime_burst_map.json`](../../reports/d1_chime_burst_map.json). Summary: 51 arc rows (48 codetection `.npy` + 3 CANFAR session dirs); all 48 `.npy` rows linked to iacobus via nickname/TNS-date alias (e.g. `johndoeii` → `johndoe_230814aaas`); **0 exact basename overlap**.

**Prior next step (superseded by map script):**

```bash
python scripts/migration/audit_arc_delta.py --stdout   # refresh arc vs iacobus counts
# Manual: map nickname ↔ TNS ↔ arc basename per burst via configs/bursts.yaml + data-manifest.csv
```

**Do not:** bulk rsync arc → iacobus without quota check (~200 G arc cap).

---

## D2 — iacobus `CHIME_canfar` archive merge

| Source | Target | Size |
|--------|--------|------|
| `~/Archives/CHIME_canfar` (iacobus) | `~/Research/CHIME_DSA_Codetections/archive/chime_canfar/` | 725 f / 2.7 G |

**Finding (2026-06-26):** zero basename overlap vs `Research/…/archive` (937 f / 178 G). Merge is **additive**, not dedupe — unique CANFAR session exports (includes 3 `analysis_*` session dirs with spaces in names).

**Completed 2026-06-27** (move-only on iacobus; no dedupe):

```bash
python scripts/migration/audit_chime_canfar.py --stdout   # pre-move inventory
# on iacobus: mv ~/Archives/CHIME_canfar ~/Research/CHIME_DSA_Codetections/archive/chime_canfar
python scripts/query_machine_inventory.py --migration-map --json | jq '.[] | select(.id=="iacobus_chime_canfar_archive")'
```

Pre-move audit: [`reports/d2_chime_canfar_inventory.csv`](../../reports/d2_chime_canfar_inventory.csv) — 725 source rows, 937 archive rows, **0 basename overlap**. Post-move verify: source path absent; target 725 f / 2.7 G.

**Inventory id:** `iacobus_chime_canfar_archive` (`status: completed`).

---

## D3 — h17 arc trash → iacobus (staging)

| Source | Target | Size |
|--------|--------|------|
| h17 `.../chime-dsa-codetections/archive/arc_trash_2026-06` (was `/data/jfaber/arc_archive_2026-06`) | iacobus `archive/arc_trash_2026-06/` | 1924 f / 36 G |

**Finding (2026-06-27 hash-map):** full sha256 on 245 `.pkl`/`.npy` (36.4 G hashed) vs iacobus `OLD_CHIME_DSA_Codetections` + `archive/chime_canfar` — **97.5% unique bytes** (20 hash duplicates / 906 M duplicate bytes). Basename overlap remains low (numeric vs nickname naming).

**Executed (2026-06-27):** rsync h17→iacobus (`03:45 UTC`); h17 source consolidated move-only into compute workspace (`07:24 UTC`).

```bash
python scripts/migration/audit_h17_arc_archive.py --stdout
bash scripts/migration/h17_to_iacobus.sh          # iacobus pull via ssh -A; --ignore-existing
python scripts/query_machine_inventory.py --migration-map --json | jq '.[] | select(.id=="h17_arc_archive_copy")'
```

Pre-move audit: [`reports/d3_h17_arc_inventory.csv`](../../reports/d3_h17_arc_inventory.csv), [`reports/d3_h17_arc_inventory.json`](../../reports/d3_h17_arc_inventory.json).

**Archive-target note (2026-06-26):** iacobus is **staging only** — the upload
target is Google Drive (`jakobtfaber@gmail.com`). This copy consolidates h17
compute artifacts onto staging ahead of D5 upload.

**Inventory id:** `h17_arc_archive_copy` (`status: completed`).

---

## D5 — iacobus → Google Drive bulk upload **COMPLETED 2026-07-13**

The upload was initially deferred on 2026-06-27, then completed directly from
iacobus on 2026-07-13. Final `rclone check --size-only` reported 0 differences
and 5,437 matching files. Google Drive holds 244.815 GiB / 5,438 objects.

| Source | Target | Size |
|--------|--------|------|
| iacobus former staging tree | `gdrive-jakob:Research/CHIME_DSA_Codetections` | 244.815 GiB |

The source tree was then moved, without deletion, to
`iacobus:~/Research/_quarantine/CHIME_DSA_Codetections-drained-20260713/`.
Full evidence:
[`HANDOFF_mbp_tailscale_ssh_iacobus.md`](HANDOFF_mbp_tailscale_ssh_iacobus.md).

**Inventory id:** `iacobus_to_gdrive` (`status: completed`).

---

## D4 — Docs: CANFAR GPU access

Local-only note through 2026-06-25; committed separately in PR (see `DATA_SOURCES.md` § CANFAR compute and GPU access).

Smoke test (verified 2026-06-25):

```bash
canfar create headless skaha/astroml-cuda:latest --gpu 1 -n gpu-smoke-test -- nvidia-smi
```

---

## Quick reference

| id | Risk if rushed | Approval needed |
|----|----------------|-----------------|
| D1 CHIME_bursts | wrong namespace / duplicate fits | yes — reconcile map first |
| D2 CHIME_canfar | none (additive, iacobus-local) | yes — move-only merge |
| D3 h17 arc trash | 36 G mostly unique (97.5% hash) | executed — iacobus staging toward gdrive |
| D4 GPU docs | none | no (docs only) |
| D5 iacobus→gdrive | completed and verified | **completed 2026-07-13** |

**Audit artifacts:** `reports/phase3_audit.json`, `reports/phase4_audit.json`, `reports/phase3_chime_basename_inventory.csv`.

**Related closeouts:** [`PHASE4_CLOSEOUT.md`](PHASE4_CLOSEOUT.md), [`PHASE5_CLOSEOUT.md`](PHASE5_CLOSEOUT.md).
