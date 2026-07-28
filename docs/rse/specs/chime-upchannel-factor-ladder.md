# CHIME upchannelization-factor ladder

**Status:** infrastructure only; authoritative generation blocked

**Factors:** 16, 32, 64, 128, 256, 512
**Scope:** CHIME/FRB input products for later autocorrelation-function analysis

No burst product or H17 job was run while adding this infrastructure.

## Gates

Authoritative generation requires, for every event:

1. a ratified dispersion-measure artifact;
2. the owner-reviewed bad-channel policy;
3. a pinned container digest and exact worker bytes;
4. a successful dry run with no output collision;
5. owner approval to start the Zach pilot.

The dispersion-measure artifact is `<dm-dir>/<event>.json`:

```json
{
  "schema": "faber2026-chime-dm-provenance-v1",
  "event": "zach",
  "instrument": "CHIME/FRB",
  "dm_pc_cm3": 262.368,
  "ratification_status": "ratified",
  "ratified_by": "owner identity",
  "ratified_at_utc": "YYYY-MM-DDTHH:MM:SSZ",
  "decision_record": "durable decision-record path or identifier"
}
```

Candidate artifacts are accepted only by `--dry-run`. Production refuses them.
There is no provisional-production switch.

## Dry run

Run inside the pinned CHIME baseband-analysis container with the staged H5 tree
mounted read-only. The command checks every input hash and H5 layout, computes
the expected frequency resolution, time resolution, shape, and output identity,
checks collisions, and writes only the planned campaign manifest.

```bash
python upchannelize_chime.py zach \
  --upchannel-factor 16 \
  --upchannel-factor 32 \
  --upchannel-factor 64 \
  --upchannel-factor 128 \
  --upchannel-factor 256 \
  --upchannel-factor 512 \
  --dm-provenance-dir /work/dm-provenance \
  --container-identity 'chimefrb/baseband-analysis@sha256:<digest>' \
  --out /work/products \
  --no-time-shift \
  --dry-run \
  --planned-manifest /work/evidence/zach-factor-ladder-plan.json
```

The output tree is immutable and separates instrument, event, dispersion
measure plus provenance hash, factor, worker hash, and container hash:

```text
instrument-chime-frb/
  event-zach/
    dm-262p368-<dm-hash>/
      u0032/
        worker-<worker-hash>/
          container-<container-hash>/
            manifest.json
            acf_input.npz
            stokes_i.npy
            frequencies_mhz.npy
            package_frequencies_mhz.npy
            source_valid.npy
```

An existing leaf directory is a collision. The worker refuses to overwrite it.

## Future execution order

1. Ratify all per-burst dispersion measures and store the artifacts.
2. Record the owner-reviewed bad-channel policy.
3. Re-run the Zach six-factor dry run.
4. Owner reviews the plan, expected shapes, storage, and exact command.
5. Generate the six Zach products.
6. Validate hashes, frequency and time resolution, shapes, nominal grid,
   source-valid mask, and product manifests.
7. Build factor-specific full-grid mappings and owner-approved effective masks.
8. Run diagnostic autocorrelation functions with explicit lag-support settings.
9. Owner reviews the Zach comparison before authorizing the 12-event,
   six-factor campaign: 72 products.

The factor ladder changes both frequency and time resolution:

| Factor | Frequency resolution | Time resolution |
|---:|---:|---:|
| 16 | 24.4140625 kHz | 0.08192 ms |
| 32 | 12.20703125 kHz | 0.16384 ms |
| 64 | 6.103515625 kHz | 0.32768 ms |
| 128 | 3.0517578125 kHz | 0.65536 ms |
| 256 | 1.52587890625 kHz | 1.31072 ms |
| 512 | 0.762939453125 kHz | 2.62144 ms |

On-pulse and off-pulse windows must be specified in physical time and converted
again for each factor. Reusing bin indices would select different time
intervals and invalidate the comparison.

## Autocorrelation-function handoff

Each completed product manifest binds the source H5, ratified dispersion
measure, factor, expected and realized grid, worker, container, command, and all
output hashes. `acf_input.npz` is directly loadable by the scintillation
pipeline.

Factor-tagged inputs must configure `analysis.upchannel_product` with the exact
manifest hash and identities. They must also use a required bad-channel mask
whose `source_valid_path` is the manifest-bound file, plus explicit
`analysis.acf.first_fit_lag`, `min_support_pairs`, and
`min_support_fraction`. The downstream cache identity includes the product
manifest bytes and the existing mask/mapping bytes. This adds no fit and changes
no scientific claim.
