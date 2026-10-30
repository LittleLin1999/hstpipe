"""Target-independent run configuration; paths are relative to the project root."""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
import json
import re


@dataclass
class BackgroundConfig:
    enabled: bool = True
    source_radius: float = 80.0
    sky_inner: float = 100.0
    spatial_widths: tuple = (25.0, 50.0, 100.0)
    continuum_bin: float = 20.0
    bootstrap_draws: int = 80
    residual_enabled: bool = True
    # Distances are offsets beyond source_radius; both strips remain INSIDE
    # sky_inner, so the primary fit has seen neither fitting nor validation sky.
    residual_fit_offsets: tuple = (4.0, 10.0)
    residual_validation_offsets: tuple = (12.0, 18.0)
    residual_detection_sigma: float = 5.0
    residual_validation_sigma: float = 3.0

    def __post_init__(self):
        self.spatial_widths = tuple(self.spatial_widths)
        self.residual_fit_offsets = tuple(self.residual_fit_offsets)
        self.residual_validation_offsets = tuple(self.residual_validation_offsets)
        if self.residual_enabled:
            a, b = self.residual_fit_offsets, self.residual_validation_offsets
            if not (len(a) == len(b) == 2 and 0 < a[0] < a[1] < b[0] < b[1]
                    and self.source_radius + b[1] < self.sky_inner):
                raise ValueError("Residual strips must be disjoint, protected, and unused by primary fit")
            if self.residual_detection_sigma < 3 or self.residual_validation_sigma < 2:
                raise ValueError("Residual detection thresholds are too permissive")
        if self.source_radius <= 0 or self.sky_inner <= self.source_radius:
            raise ValueError("sky_inner must exceed the protected source radius")
        if not self.spatial_widths or min(self.spatial_widths) <= 0:
            raise ValueError("Spatial widths must be positive")
        if self.continuum_bin <= 0 or self.bootstrap_draws < 10:
            raise ValueError("Invalid background sampling")


@dataclass
class ExtractionConfig:
    """Separate boxcar, optimal, and training heights shared by all branches."""

    enabled: bool = True
    height_native_pixels: float = 11.0
    # G140L PCTAB contains 140 pixels; 61 would silently select 21 in HSTCAL.
    optimal_height_native_pixels: float = 140.0
    profile_height_native_pixels: float = 151.0
    optimal_min_profile_fraction: float = 0.8
    optimal_min_information_fraction: float = 0.5
    center_method: str = "continuum"
    center_arcsec: float = 0.0
    trace_tolerance_pixels: float = 0.35
    # Legacy x1d diagnostic geometry; retained for old configs, unused with OMIT.
    residual_offset_native_pixels: float = 112.5
    residual_size_native_pixels: float = 25.0
    flux_mode: str = "point"
    source_model: str | None = None

    def __post_init__(self):
        import math

        if not all(
            math.isfinite(v)
            for v in (
                self.height_native_pixels,
                self.optimal_height_native_pixels,
                self.profile_height_native_pixels,
                self.optimal_min_profile_fraction,
                self.optimal_min_information_fraction,
                self.center_arcsec,
                self.trace_tolerance_pixels,
                self.residual_offset_native_pixels,
                self.residual_size_native_pixels,
            )
        ):
            raise ValueError("Extraction options must be finite")
        if not (0 < self.optimal_min_profile_fraction <= 1
                and 0 < self.optimal_min_information_fraction <= 1):
            raise ValueError("Optimal coverage thresholds must be in (0, 1]")
        if self.center_method not in ("header", "continuum"):
            raise ValueError("center_method must be header or continuum")
        if not (self.height_native_pixels <= self.optimal_height_native_pixels
                <= self.profile_height_native_pixels):
            raise ValueError("Require boxcar <= optimal <= profile fitting height")
        if self.height_native_pixels < 2 or self.trace_tolerance_pixels <= 0:
            raise ValueError("Extraction height must be at least two native pixels")
        if self.residual_size_native_pixels <= 0 or self.residual_offset_native_pixels <= 0:
            raise ValueError("Residual background windows must have positive geometry")
        if self.flux_mode not in ("point", "diffuse", "morph"):
            raise ValueError("flux_mode must be point, diffuse, or morph")
        if self.flux_mode == "morph" and not self.source_model:
            raise ValueError("morph flux calibration requires a source_model")


@dataclass
class Config:
    project_dir: str = "."
    input_glob: str = "raw/*_flt.fits"
    mast_dir: str = "raw"
    reference_dir: str = "calibration_refs"
    control_dir: str = "intermediate/flt/registration/official"
    output_dir: str = "reduction"
    software_dir: str | None = None
    name: str | None = None
    spatial_scale: float = 0.0123
    spatial_limits: tuple | None = None
    pixfrac: float = 1.0
    source_offsets: dict = field(default_factory=dict)
    prepare_controls: bool = False
    background: BackgroundConfig = field(default_factory=BackgroundConfig)
    mast_background: BackgroundConfig = field(default_factory=BackgroundConfig)
    extraction: ExtractionConfig = field(default_factory=ExtractionConfig)

    def __post_init__(self):
        if isinstance(self.background, dict):
            self.background = BackgroundConfig(**self.background)
        if isinstance(self.mast_background, dict):
            self.mast_background = BackgroundConfig(**self.mast_background)
        if isinstance(self.extraction, dict):
            self.extraction = ExtractionConfig(**self.extraction)
        if self.name is not None and not re.fullmatch(r"[A-Za-z0-9_.+-]+", self.name):
            raise ValueError("name must be a filename-safe label")
        if self.spatial_scale <= 0 or not 0 < self.pixfrac <= 1:
            raise ValueError("Invalid output sampling")
        if self.spatial_limits is not None:
            if (
                len(self.spatial_limits) != 2
                or self.spatial_limits[0] >= self.spatial_limits[1]
            ):
                raise ValueError("Invalid spatial_limits")

    def path(self, value):
        p = Path(value).expanduser()
        return (
            p.resolve() if p.is_absolute() else (Path(self.project_dir) / p).resolve()
        )

    @property
    def output(self):
        return self.path(self.output_dir)

    @classmethod
    def read(cls, filename):
        filename = Path(filename).resolve()
        values = json.loads(filename.read_text())
        root = Path(values.get("project_dir", ".")).expanduser()
        values["project_dir"] = str(
            (filename.parent / root).resolve() if not root.is_absolute() else root
        )
        return cls(**values)

    def serializable(self):
        return asdict(self)

    def resolve_name(self, exposures):
        """Derive the product prefix from the science FITS target and grating.

        An explicit name is accepted only when it agrees with the headers, so
        a copied config cannot silently label another target's products.
        """
        identities = {
            (str(e.p.get("TARGNAME", "")).strip(), str(e.p.get("OPT_ELEM", "")).strip())
            for e in exposures
        }
        if len(identities) != 1 or not exposures:
            raise ValueError("Inputs have inconsistent target or OPT_ELEM headers")
        target, opt_elem = identities.pop()
        prefix = f"{target}_{opt_elem}"
        if not target or not opt_elem or not re.fullmatch(r"[A-Za-z0-9_.+-]+", prefix):
            raise ValueError("FITS target or OPT_ELEM is unsuitable for a product name")
        if self.name is not None and self.name != prefix:
            raise ValueError(f"Configured name {self.name!r} differs from FITS {prefix!r}")
        self.name = prefix
        return prefix

    def product_name(self, kind):
        """One source of truth for the three final science-product names."""
        if self.name is None:
            raise ValueError("Resolve the FITS target and OPT_ELEM before naming products")
        suffixes = {
            "s2d_drz": "_s2d_drz.fits",
            "x2d_coadd": "_x2d_coadd.fits",
            "x1d": "_x1d.fits",
        }
        return self.name + suffixes[kind]
