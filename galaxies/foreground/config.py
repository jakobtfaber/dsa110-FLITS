"""Configuration and target list for foreground galaxy search."""

from astropy.cosmology import Planck18

# Cosmology to use throughout the module
COSMO = Planck18

# Default search parameters
DEFAULT_IMPACT_KPC = 100.0
DEFAULT_CLUSTER_IMPACT_KPC = (
    5000.0  # 5 Mpc: clusters intercept sightlines far beyond galaxy-scale impacts
)
DEFAULT_Z_EPS = 0.01  # Redshift buffer for foreground search
MIN_Z_SEARCH = 0.005  # Minimum redshift cutoff to avoid infinite search cone
# Minimum credible photometric redshift: DESI Legacy photo-z (VII/292/north) pin
# unreliable sources near z~0 (z=0.001-0.008 floor noise). Applies only to rows
# carrying a photo-z error; spec-z catalogs (NED, GLADE+) are exempt so genuine
# nearby galaxies (e.g. UGC 06371, z=0.009) survive.
FOREGROUND_PHOTOZ_FLOOR = 0.01
# Velocity-offset window (km/s) below which a spec-z neighbour is host/local
# ambiguous, not a clean foreground. dv = c*(z_frb - z_gal)/(1 + z_frb); a galaxy
# inside this window sits at the host's recession velocity (group member / host
# peculiar velocity ~200-500 km/s) and cannot be cleanly attributed as an
# intervening system. Set at the upper edge of typical group dispersions.
FOREGROUND_AMBIGUITY_KMS = 500.0
SPEED_OF_LIGHT_KMS = 299792.458
# Spec-z catalogs: redshift is spectroscopic / distance-derived, so any z_err they
# carry (e.g. GLADE+'s generic 0.015 floor) is NOT a per-object photo-z and must
# not gate via the photo-z error path. Matched on the `catalog` column substring;
# DESI VII/292/north is the lone photo-z catalog.
SPEC_Z_CATALOG_SUBSTRINGS = ("ned", "vii/291", "glade", "sdss", "desi_dr1")
PHOTO_Z_CATALOG_SUBSTRINGS = ("vii/292", "legacy_dr9_photoz")
# The 12 FRB sightlines in our sample. ``None`` means that no host redshift is
# established; acquisition may return candidates, but must not call them
# foreground/background until a redshift posterior is supplied.
# Format: (name, RA, Dec, z_frb)
TARGETS: list[tuple[str, str, str, float | None]] = [
    ("Zach", "20h40m47.886s", "+72d52m56.378s", 0.0430),
    ("Whitney", "08h58m52.92s", "+73d29m27.0s", 0.4790),
    ("Oran", "21h12m10.760s", "+72d49m38.20s", 0.3005),
    ("Isha", "04h45m38.64s", "+70d18m26.6s", 0.2505),
    ("Wilhelm", "21h00m31.09s", "+72d02m15.22s", None),
    ("Phineas", "11h51m07.52s", "+71d41m44.3s", 0.2710),
    ("Freya", "05h52m45.12s", "+74d12m01.7s", None),
    ("Hamilton", "20h20m08.92s", "+70d47m33.96s", 0.3024),
    ("Mahi", "02h39m03.96s", "+71d01m04.3s", None),
    ("Chromatica", "20h50m28.59s", "+73d54m00.0s", 0.0740),
    ("Casey", "11h19m56.05s", "+70d40m34.4s", 0.2870),
    ("Johndoeii", "22h23m53.94s", "+73d01m33.26s", 0.5535),
]

# Diagnostic DM-z distributions from analysis/scripts/dm_redshift_inference.py.
# These are not established host redshifts and must never enter a point-estimate
# DM budget. They permit probability-labeled candidate triage only.
TARGET_DM_REDSHIFT_ESTIMATES = {
    "Wilhelm": {"z16": 0.4102, "z50": 0.5492, "z84": 0.7048},
    "Freya": {"z16": 0.7256, "z50": 0.9173, "z84": 1.1441},
    "Mahi": {"z16": 0.7365, "z50": 0.9317, "z84": 1.1631},
}

# Catalog identifiers for Vizier
VIZIER_CATALOGS = {
    "GLADE+": "VII/291/gladep",  # GLADE+ (2022); Vizier renamed table glade -> gladep
    "DESI_DR8_NORTH": "VII/292/north",  # High-Dec northern sky coverage
    "SDSS_DR12": "V/147/sdss12",  # Stable spectroscopic/photometric catalog
    "GSC242": "I/353/gsc242",  # Guide Star Catalog 2.4.2 (star/galaxy morphology classification)
    "CATWISE2020": "II/365/catwise",  # CatWISE2020 IR proper motion & photometry
    "UNWISE": "II/363/unwise",  # unWISE deep co-adds
}

# Opt-in extra search engines (TAP-backed). Disabled by default so run_search()
# output is byte-for-byte unchanged unless explicitly enabled.
EXTRA_SEARCH_ENGINES = {"DESI_DR1": "desi_dr1.zpix"}
LEGACY_DR9_PHOTOZ_ROOT_URL = "https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr9"
LEGACY_DR9_PHOTOZ_CACHE_ENV = "FLITS_LEGACY_DR9_SWEEP_CACHE"

# Opt-in enricher catalogs (cross-matched photometry, not foreground search).
ENRICHER_CATALOGS = {
    "DESI_LS_DR10": "ls_dr10.tractor",
    "ALLWISE": "II/328/allwise",
    "CATWISE2020": "II/365/catwise",
    "UNWISE": "II/363/unwise",
    "GSC242": "I/353/gsc242",
    "GALEX_AIS": "II/335/galex_ais",
    "2MASS_XSC": "VII/233/xsc",
}

ENABLE_EXTRA_ENGINES = False
ENABLE_LEGACY_DR9_PHOTOZ = False
ENABLE_ENRICHERS = False

# All-sky galaxy-cluster catalogs. Only all-sky catalogs reach the sample's
# declination (+70..+74); SDSS-based WHL/redMaPPer do not (see
# docs/rse/specs/research-foreground-galaxies-sightlines.md). PSZ2 reports M500 as
# MSZ (1e14 Msun); MCXC/MCXC-II report M500 (1e14) + R500 (Mpc).
CLUSTER_VIZIER_CATALOGS = {
    "WEN_HAN_2024": "J/ApJS/272/39/table2",
    "PSZ2": "J/A+A/594/A27/psz2",
    "MCXC": "J/A+A/534/A109/mcxc",
    "MCXC_II": "J/A+A/688/A187/mcxcii",  # MCXC-II (A&A 688, A187; arXiv:2402.01538), live-confirmed
}

# Pinned acquisition contract. Depth values are approximate catalog selection
# limits, not local exposure-map completeness; per-position depth remains a
# required downstream qualification field.
SURVEY_CONTRACT = {
    "NED": {"release": "live TAP snapshot", "depth": "heterogeneous", "role": "discovery"},
    "GLADE+": {"release": "VII/291/gladep", "depth": "heterogeneous", "role": "discovery"},
    "DESI_DR8_NORTH": {
        "release": "VII/292/north",
        "depth": "Legacy DR8 North source selection",
        "role": "photo-z discovery",
    },
    "SDSS_DR12": {"release": "V/147/sdss12", "depth": "DR12 selection", "role": "discovery"},
    "GSC242": {"release": "I/353/gsc242", "depth": "catalog selection", "role": "classification"},
    "CATWISE2020": {
        "release": "II/365/catwise",
        "depth": "CatWISE2020 selection",
        "role": "classification",
    },
    "UNWISE": {"release": "II/363/unwise", "depth": "unWISE selection", "role": "classification"},
    "CLUSTERS": {
        "release": "Wen & Han 2024 + PSZ2 + MCXC + MCXC-II",
        "depth": "Wen-Han M500 >= 4.7e13 Msun; other catalog-specific selections",
        "role": "cluster discovery",
    },
}
# Keep a cluster when impact <= this multiple of its own r200 (research: meaningful
# cluster DM needs the sightline within ~1-2 r200).
CLUSTER_R200_FACTOR = 2.0
ENABLE_CLUSTER_ENGINE = True
