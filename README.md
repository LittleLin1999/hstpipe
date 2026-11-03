# hstpipe

Reduce calibrated Hubble Space Telescope observations into STIS spectra or
WFC3/UVIS images. Original exposures stay unchanged. Each reduction keeps its
configuration, diagnostics, and processing record alongside the products.

| Your task | Start here |
| --- | --- |
| Install the package | [Installation](docs/INSTALLATION.md) |
| Reduce UVIS images | [UVIS user guide](docs/uvis.md) |
| Understand a UVIS report or approve a reviewed image | [Review and publication](docs/uvis.md#review-and-publication) |
| Reduce STIS spectra | [STIS workflow](docs/WORKFLOWS.md#stis-reduction) |
| Remove regenerable UVIS diagnostics | [Cleanup](docs/uvis.md#cleanup) |
| Understand the code | [Architecture](docs/ARCHITECTURE.md) |

## Install

The full setup is tested on **Linux with Python 3.11**. It requires git, CMake,
a C compiler, pkg-config, and CFITSIO development libraries. See the
[installation guide](docs/INSTALLATION.md) for prerequisites and instrument-only
installation.

```bash
git clone https://github.com/LittleLin1999/hstpipe.git
cd hstpipe
bash setup.sh "$HOME/.venvs/hstpipe"
source "$HOME/.venvs/hstpipe/bin/activate"
```

The setup installs both instrument dependencies and the patched HSTCAL needed
for STIS. It does not download observations or calibration references. Run
`bash setup.sh --check` to check prerequisites without installing anything.
The complete native setup has not been validated on macOS or Windows.

## Keep observations separate from code

Run the installed package against a data directory of your choice:

```text
observations/
├── raw/
│   ├── stis/              # Calibrated FLT and matching X2D files
│   └── uvis/              # Calibrated, CTE-corrected FLC files
├── calibration/
│   ├── stis/
│   └── uvis/
├── configs/
└── redux/                 # Reduction outputs
```

Place inputs directly in each raw directory, retaining their MAST filenames.
Calibration references must match the observations. UVIS and STIS automatically
retrieve missing header-named files through HST CRDS and reuse local copies. See
[reference locations and offline settings](docs/INSTALLATION.md#download-calibration-references). Paths in YAML are relative to that YAML file. Command-line paths are
relative to the current directory.

## Reduce UVIS images

From the observations directory, with the environment activated:

```bash
export iref=/absolute/path/to/observations/calibration/uvis/
hstpipe uvis inspect --input raw/uvis
hstpipe uvis reduce --input raw/uvis --output redux --dry-run
hstpipe uvis reduce --input raw/uvis --output redux --resume
```

Add `--target TARGET --filter F336W` to select one group. A dry run prints the
plan; it does not validate the scientific result.

The pipeline measures relative alignment using sources on every available CCD,
then makes a provisional mosaic. Gaia/Legacy absolute alignment is propagated
back to working FLCs through TweakBack before the final AstroDrizzle pass.
Final defaults are **0.04 arcsec/pixel** and **pixfrac 0.8**.

Open `redux/{target}/uvis/products/{target}_uvis_report.html` first. It works
offline and is smaller than 1,000,000 bytes. Published images are named
`{target}_{filter}_drc.fits`; filters share one target directory.

A completed command does not guarantee science readiness. A blocked reduction
keeps diagnostics. **Review image QA** requires inspection; eligible field-source
flux exceptions can be explicitly approved without changing image pixels.
See [what to review and what happens after approval](docs/uvis.md#review-and-publication).

## Reduce STIS spectra

The supported mode is **STIS FUV-MAMA/G140L**. It uses calibrated FLTs, matching
X2Ds, reference tables, and the supplied patched HSTCAL build.

```bash
hstpipe stis inspect --input raw/stis
hstpipe stis reduce --config configs/stis.yaml --dry-run
hstpipe stis reduce --config configs/stis.yaml --resume
```

See the [STIS workflow](docs/WORKFLOWS.md#stis-reduction) for a complete YAML
example, control-image requirements, and product definitions. Outputs remain
under `redux/{target}/stis/{grating}/`. Other STIS modes, ACS, and COS are not
currently supported.

## Use from Python

```python
from hstpipe.uvis import reduce_uvis

plan = reduce_uvis("configs/uvis.yaml", dry_run=True)
results = reduce_uvis("configs/uvis.yaml")
```

`hstpipe.stis.reduce_stis` has the same configuration interface. Both accept a
YAML filename, a `ReductionConfig`, or an equivalent dictionary.

## More detail

- [Workflows, output definitions, and restart behavior](docs/WORKFLOWS.md)
- [UVIS source detection and multi-chip alignment](docs/uvis-full-field.md)
- [UVIS relative-fit models and acceptance](docs/uvis-relative.md)
- [DQ, sky, CR, and photometric policies](docs/uvis-handbook.md)
- [Historical validation records](docs/VALIDATION.md)

Historical test counts and target results describe particular runs. For current
science readiness, use the target's report and `manifest.json`.
