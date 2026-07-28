# ACF RFI provenance and lag-support qualification

**Date:** 2026-07-28
**Scope:** software qualification only

The ACF preparation paths now accept an explicit full-grid-to-compact-row
mapping. The mapping uses integer row identities; frequency values only verify
the declared ordering. Input, full and compact axes, source-valid rows,
owner-map bytes, effective-mask bytes, mapping bytes, event, and instrument
must all agree. A configured required mask fails closed.

The effective mask is applied before frequency-grid regularization and
frequency averaging. Bursts without configured mask provenance keep legacy
unmasked behavior.

Both ACF implementations expose surviving pair counts and the fraction of
possible pairs surviving at each lag. `analysis.acf.min_support_pairs` defaults
to `2`; `analysis.acf.min_support_fraction` defaults to `0.0`, preserving the
existing numerical floor. Declaring stricter limits removes unsupported lag
points. These defaults are software behavior, not scientific admission.

Cache names bind current source-valid, owner-map, effective-mask, mapping, and
provenance bytes. Cached spectrum and ACF loads revalidate the configured
artifacts.

Tests cover candidate upchannelization factors 16, 32, 64, 128, 256, and 512;
frequency averaging; partial and fully masked blocks; extreme-value
non-leakage; exact lag values, pair counts, and weights; missing or altered
provenance; axis/order failures; cache invalidation; and both preparation
entry points.

No burst was replayed. No result, mask, fit, or manuscript claim was promoted.
Scientific use still requires a provenance-valid mapping for the exact input,
an owner-approved effective mask, a declared support threshold, replay, and
owner review.
