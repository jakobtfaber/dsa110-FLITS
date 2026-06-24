# Visual component census (Jakob, 2026-06-23)

Per-band temporal sub-pulse counts, marked by eye on `component_screen/*_components.png`
(dedispersed waterfall + frequency-collapsed profile, on-pulse cropped, no model
overlay) via the interactive `component_deck.html`. Annotated PNGs in
`component_screen/annotated/`. These visual counts are the INPUT; the
marginalize-gain multi-component evidence ladder (lnZ) confirms them numerically.

| burst | CHIME | DSA | class | notes |
|---|---|---|---|---|
| hamilton | 5 | 1 | multi-CHIME | C1-C5 on a narrow spike cluster ~11.5-13 ms; DSA single narrow spike |
| zach | 2 | 3-4 | multi both | C1+C2; DSA D1 main + structured secondary (D2/D3/D4 ~11-13 ms) |
| phineas | 3 | 3 | multi both | C1-C3 across the broad CHIME envelope; D1-D3 |
| johndoeII | 2 | 2 | multi both | C1+C2 (tail shoulder); D1 main + D2 ~1.4 ms |
| whitney | 2 | 1 | multi-CHIME | C1+C2 ~1.5-2.5 ms; DSA single narrow spike |
| isha | 1 | 1 | single | DSA shoulder read as scattering tail, not 2nd pulse |
| mahi | 1 | 1 | single | DSA shoulder = tail |
| oran | 1 | 1 | single | DSA shoulder = tail |
| casey | 1 | 1 | single | clean narrow singles both bands |
| chromatica | 1 | 1 | single | clean |
| freya | 1 | 1 | single | clean (control) |
| wilhelm | 1 | 1 | single | DSA broad single |

**5 multi-component bursts** (hamilton, zach, phineas, johndoeII, whitney);
**7 single-component** (isha, mahi, oran, casey, chromatica, freya, wilhelm).

The key disambiguation the visual screen settled: isha/mahi/oran show a DSA-band
shoulder that the chi2-tension alone could not distinguish from a second pulse;
marked as 1 each (tail, not component).

A naive arrow-color detector under-counted tightly-clustered marks (e.g. zach DSA
secondary, hamilton CHIME cluster) -- counts above are the by-eye read of the
annotated PNGs, which supersede the detector.
