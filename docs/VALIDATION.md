# Code-only validation

2026-10-30: editable installation succeeded on erebus with Python 3.11.
62 tests pass (53 preserved STIS tests and 9 infrastructure/adapter tests).
Five existing low-sample background variance warnings remain.
The test process used the existing Drizzle 3 source installation because the
server's default environment has an older incompatible Drizzle API.

Tests covered masked extraction, background algorithms, coaddition math,
calibration conversions, hstpipe public imports, public imports, installed CLI,
YAML paths, raw/redux separation, product names, file inventory, duplicate
contents, subarray metadata, cache invalidation, and mocked DrizzlePac dispatch.
The recorded AST comparison confirms only the documented path literals changed
in the preserved STIS modules. This is not a real-data numerical regression.

No observational reduction was run. No real UVIS alignment, cosmic-ray rejection,
sky subtraction or photometric validation was performed. Package-level imports
were checked without loading DrizzlePac. See ARCHITECTURE.md for remaining
stage-cache and directory migration work.

Compatibility namespace removal: all STIS tests now import hstpipe.stis;
only the hstpipe console entry point is installed. No science module changed.
