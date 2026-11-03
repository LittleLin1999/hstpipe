# Reduction workflows

[Home](../README.md) · [UVIS guide](uvis.md) · [Installation](INSTALLATION.md)

Use this page for the STIS configuration and shared restart behavior. UVIS has a
[separate user guide](uvis.md) covering its report, image review, and cleanup.

## Choose an instrument

| Instrument | Input | Output directory |
| --- | --- | --- |
| STIS FUV-MAMA/G140L | Calibrated FLTs, matching X2Ds, reference tables | `redux/{target}/stis/{grating}/` |
| WFC3/UVIS | Calibrated, CTE-corrected FLCs and distortion references | `redux/{target}/uvis/` |

## STIS reduction

First inspect the exposures. The command prints a JSON inventory and compatible
target/grating groups:

```bash
hstpipe stis inspect --input /path/to/observations/raw/stis
```

Save the following as `observations/configs/stis.yaml`. Replace `J1022+0841` with
the target name in your FITS headers:

```yaml
input: ../raw/stis
output: ../redux
instrument: stis
target: J1022+0841
element: G140L
options:
  reference_dir: ../calibration/stis
  prepare_controls: true
  extraction:
    height_native_pixels: 11
    optimal_height_native_pixels: 140
    profile_height_native_pixels: 151
    flux_mode: point
```

`prepare_controls: true` generates the official x2d controls using the patched
HSTCAL installed by the setup script. To reuse existing controls, set
`options.control_dir` to their directory and `options.prepare_controls: false`.
That directory must contain `*_flt_omit.fits`, `*_flt_perform.fits`, and
`*_flt_mapping.f64` for each input FLT. Use reference tables that match the
observations; changing calibration versions changes the reduction.

Preview the selected inputs and output location, then start processing:

```bash
hstpipe stis reduce --config /path/to/observations/configs/stis.yaml --dry-run
hstpipe stis reduce --config /path/to/observations/configs/stis.yaml
```

`--dry-run` checks input selection and prints the plan; it does not execute science
stages or establish that all calibration dependencies are ready.

To continue a run or deliberately start it again:

```bash
# Reuse verified outputs; rerun stages whose inputs or dependencies changed.
hstpipe stis reduce --config /path/to/observations/configs/stis.yaml --resume

# Preserve the existing run in a timestamped backup and rebuild.
hstpipe stis reduce --config /path/to/observations/configs/stis.yaml --force
```

YAML paths are relative to the YAML file. Explicit command-line paths are relative
to your current directory. Separate incompatible instrument setups into separate
input selections; `--setup` adds an output label and does not filter exposures.

### STIS products

For the example above, outputs are written under
`redux/J1022+0841/stis/G140L/`:

| Location | Contents |
| --- | --- |
| `products/*_s2d_drz.fits` | Combined spectrum from the FLT/drizzle branch |
| `products/*_x2d_coadd.fits` | Combined spectrum from the x2d branch |
| `products/*_x1d.fits` | FLT/DRZ/X2D tables containing boxcar and optimal extractions |
| `diagnostics/` | Background, combination, extraction, and flux-calibration checks |
| `intermediate/` | Per-exposure models and intermediate products |
| `logs/` | Reduction and stage logs |
| `config.yaml`, `manifest.json` | Effective configuration, input identities, and provenance |

Inspect diagnostics before interpreting the spectra. For optimal extractions,
use `DQ_OPT`; `DQ` is the boxcar mask (`DQ_BOX`). See
[stage APIs](#python-stage-apis) for the complete product layout.

## UVIS reduction

Follow the [UVIS quick start](uvis.md#run-a-reduction). All available SCI extensions
are detected automatically and paired with their ERR, DQ, and distortion WCS.
Gaia/Legacy catalogs can help identify sources, but relative fitting uses measured
HST centroids. Absolute calibration follows provisional drizzle; TweakBack then
updates the working FLCs before final AstroDrizzle.

Filters share a target directory. `intermediate/` has only `catalogs`,
`background`, `alignment`, `drizzle`, and `working` stages as needed. Diagnostics
and logs are flat. See the [output tree](uvis.md#output-layout).

## Resume, rebuild, or clean up

| Action | Command option | Effect |
| --- | --- | --- |
| Preview | `reduce --dry-run` | Show selected inputs and paths; no science processing. |
| Continue | `reduce --resume` | Reuse stages whose dependencies and output identities still match. |
| Rebuild | `reduce --force` | Preserve the previous run in a timestamped recovery location and rebuild. |
| Prune UVIS comparison images | `uvis cleanup --output redux --apply` | Delete eligible regenerable comparisons; retain final and review images. |

Changes to inputs, scientific parameters, code, dependencies, or native
executables can invalidate affected caches. A report edit is not scientific
revalidation. Cleanup does not rewrite receipts to pretend missing files exist;
a later resume may regenerate removed diagnostics. See [cleanup scope](uvis.md#cleanup).

## Python stage APIs

`hstpipe.stis.reduce_stis` and `hstpipe.uvis.reduce_uvis` accept a YAML filename,
a `ReductionConfig`, or a dictionary. The CLI provides the same top-level
`inspect` and `reduce` operations.

`hstpipe.stis.pipeline.stages()` exposes `prepare`, `background`, `joint_drizzle`,
`mast_coadd`, `validate`, and `extract`. These retain their native configuration
and input models. See the function docstrings for required arguments.

STIS keeps its two independent 2D branches and consolidated FLT/DRZ/X2D extraction
tables. The X1D `DQ` column aliases `DQ_BOX`; use `DQ_OPT` for optimal extraction.
No UVIS maintenance or report operation changes STIS scientific processing.
