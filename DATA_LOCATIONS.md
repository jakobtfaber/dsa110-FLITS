# Data Locations for CHIME-DSA Co-Detection Project

**Raw-data authority (2026-07-24):** **h17** (`lxd110h17`) — `/data/Faber2026/data/`. Every raw CHIME/FRB baseband and DSA-110 filterbank used by this project is pulled from h17 and from nowhere else. Full inventory with checksums: Faber2026 `docs/rse/ops/raw-data-provenance.md`; ledger at `h17:/data/Faber2026/provenance/h17-source-data-migration-20260721.json`.

**Host roles:** jakob-mbp (code), **h17** (raw-data authority + compute), and arc/CANFAR (CHIME baseband upstream origin; `.npy` replica). Retired hosts h23, hpcc, and dsacamera are read-only quarantine references only.

Inventory: [`machine_inventory.yaml`](machine_inventory.yaml) · Query: [`scripts/query_machine_inventory.py`](scripts/query_machine_inventory.py)

## Code (GitHub canonical)

| Host | Path |
|------|------|
| **jakob-mbp** | `~/Developer/repos/github.com/jakobtfaber/dsa110-FLITS` |
| **GitHub** | https://github.com/jakobtfaber/dsa110-FLITS |

Do not develop on hpcc, arc checkout, or h23 trees. h17 may hold an optional clone for docker workflows.

## Raw data (h17 authority)

Raw input for all twelve co-detected bursts:

| Instrument | Path on h17 |
|---|---|
| CHIME/FRB singlebeam baseband | `/data/Faber2026/data/chime-frb/<burst>/singlebeam_<chime_event_id>.h5` |
| DSA-110 filterbank, Stokes I | `/data/Faber2026/data/dsa-110/<burst>/<dsa_observation_id>_dev_polcal_I.fil` |

`/data` is a separate 13 TB volume on h17, not part of its root filesystem — a search rooted at `/` with `-xdev` misses all of it.

Upstream origins, recorded for auditing and **not** access paths: CANFAR `arc:projects/chime_frb/data/chime/baseband/processed/…` for CHIME baseband; `dsa-storage` (`dsa-storage.ovro.pvt`) and `h23` for DSA filterbanks.

Stokes Q, U and V were not migrated — h17 is authoritative for total intensity only. Polarisation work needs the four-Stokes sets staged there first.

Derived `.npy` waterfalls under `h17:/home/ubuntu/flits-runs/data/` are **not** raw data: the dispersion measure is baked into the filename and the array, so they cannot be used to revalidate a dispersion measure.

## Burst `.npy` for fits (arc + local replica)

| Host | Path | Role |
|------|------|------|
| **arc** | `arc:home/jfaber/baseband_morphologies/chime_dsa_codetections/data/DSA_bursts` | CANFAR primary (~2.9G) |
| **arc** | `.../CHIME_bursts` | Separate CHIME/FRB namespace |
| **jakob-mbp** | `~/Data/Faber2026/chimefrb/CHIME_bursts` | Current CHIME/FRB full-resolution replica |
| **jakob-mbp** | `~/Data/Faber2026/dsa110/DSA_bursts` | Current DSA-110 full-resolution replica (moved from `~/Developer/dsa110-local-data/DSA_bursts` 2026-06-30) |

Current dual-band producers require both local roots explicitly. They do not
search for CHIME/FRB products under the DSA-110 root.

## Compute workspace (h17)

Distinct from the raw-data root above. `/data/Faber2026/data/` is authoritative raw input; the paths below are working space.

| Path | Size | Notes |
|------|------|-------|
| `/data/Faber2026/` | ~24G | **Raw-data authority** — see the raw data section above; also holds `evidence/` and `provenance/` |
| `/data/research/astrophysics/frbs/chime-dsa-codetections` | ~72G | Compute + artifact cache (not a git source of truth). Held the raw files before the 2026-07-21 migration; paths into it for raw input are stale |
| `.../dsa110-FLITS/` | clone | Canonical FLITS checkout on h17 (`jakobtfaber/dsa110-FLITS`) |
| `.../upchan_codetections` | ~1.8G | Upchan products (12-target table in baseband_recovery worker) |
| `.../archive/arc_trash_2026-06` | 36G | Cold arc trash copy |
| `.../scripts/h17_codetections/` (in clone) | — | Promoted compute workers (see that README) |

Full layout, promotion rules, and sync commands:
[`docs/infrastructure/H17_WORKSPACE.md`](docs/infrastructure/H17_WORKSPACE.md).

`/data/jfaber/` is empty of codetection products as of 2026-06-27.

## Results library (processed fit products)

**Canonical local inventory (jakob-mbp):** `~/Data/Faber2026/results-library/`  
Override: `FABER2026_RESULTS_LIBRARY`.

| Role | Path |
|------|------|
| Library root | `~/Data/Faber2026/results-library/` |
| Index | `$FABER2026_RESULTS_LIBRARY/INDEX.md` |
| Catalog (git) | sibling analysis repo `../analysis/scripts/results_library_catalog.yaml` |
| Materialize | `python3 ../analysis/scripts/materialize_results_library.py` |

Phase B moves bulk campaign `results/`, `_a1_fits/`, `joint_json/`, and `pipeline/results` into the library; in-repo paths are local symlinks (gitignored). Driver `.py` stays in FLITS / parent `analysis/`. Overleaf `figures/` stay in Faber2026 git (link-only).

Full posteriors / dynesty samples remain under `$FLITS_RUNS` (see compute-scratch), not the results library.

## Retired hosts (quarantine / read-only)

Move-only policy; restore one-liners in each host's `_quarantine/README.md`.

| Host | Quarantine path | Status |
|------|-----------------|--------|
| **h23** | `/media/ubuntu/ssd/_quarantine/jfaber-drain-20260625/` | Partial: archive, burstprop_paper, chime_dsa_codetections (~137G). Residual at `jfaber/`: nihari, tools, dsa110-continuum, frb_inventory, scratch |
| **h23** | `/dataz/dsa110/T3/` | **Not quarantined** — 59T raw pipeline; leave on source |
| **hpcc** | `/home/jfaber/_quarantine/flits-20260625` | Full flits tree quarantined 2026-06-25; JSON artifacts on jakob-mbp |
| **dsacamera** | — | Decommissioned; no codetection content |
## CANFAR arc compute

Storage: `vos`/`vls` with `~/.ssl/cadcproxy.pem`. Compute: `canfar` CLI (`canfar create --gpu N`). Notebook sessions without `--gpu` are CPU-only. Live A100 MIG smoke test passed 2026-06-25.

## Related repos

- subhalos: https://github.com/jakobtfaber/subhalos — archived 2026-07-05; consolidated into frb-foreground-halos (June 2026), itself integrated into `galaxies/foreground/vo/`
- frb-foreground-halos: https://github.com/jakobtfaber/frb-foreground-halos — archived 2026-07-05; integrated into `galaxies/foreground/vo/` (PR #123). Physical results data: `~/Data/frb-foreground-halos/results/` (symlink pattern per the ~/Data convention)
- los_halos: https://github.com/jakobtfaber/los_halos — archived 2026-07-05; VO-TAP pipeline integrated into `galaxies/foreground/vo/`
