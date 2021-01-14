## AcidGenomes 0.2.0 (UNRELEASED)

### New functions

- Added new genome download functions, migrated from previous Python approach
  defined in py-koopa package.
- Switched from cli to AcidCLI for interactive messages.

### Major changes

- In metadata columns, renamed "ID" to "Id", for stricter lower camel case
  formatting.
- Now ensuring GenomicRanges gets attached as a dependency.
- Now using BiocFileCache (via `pipette::cacheURL` internally) to automatically
  cache GFF/GTF files when used in `makeGRangesFromGFF`.
- Renamed `ignoreTxVersion` to simply `ignoreVersion`, where applicable.
  We want this setting to also apply at gene level.
- No longer using `forceDetach` internally to force unload some Bioconductor
  annotation packages.
- `GRanges` `mcols` now return with `tx` prefix instead of `transcript`.
- `GRanges` `mcols` now use strict camel case formatting
  (e.g. `geneId` instead of `geneID`).

## AcidGenomes 0.1.1 (2020-10-12)

### Minor changes

- Relaxed stringency of internal organism name checks. Applies to
  `makeGene2SymbolFromEnsembl` for example, which is causing running examples
  to fail in pointillism without a fix.
- Updated minimum dependency versions.

## AcidGenomes 0.1.0 (2020-10-07)

Initial release, consisting of functions migrated from basejump.
