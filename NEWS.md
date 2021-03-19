## AcidGenomes 0.2.10 (2021-03-19)

### Minor changes

- `HGNC` now returns columns with split values as `CharacterList`, instead of
  as character strings containing "|".

## AcidGenomes 0.2.10 (2021-03-15)

### Minor changes

- `mapHumanOrthologs`: Hardened mouse-to-human matching.
- `makeGRangesFromEnsembl`: No longer hard-coding minimum release version check
  at 87, in case older releases are ported to AnnotationHub in a future release.
- Revert back to `ignoreVersion = TRUE` by default for genome annotation
  importers, as this is typically what users expect by default.

## AcidGenomes 0.2.9 (2021-03-03)

### Minor changes

- `Gene2Symbol` functions now preserve metadata, as expected. This was
  causing pointillism package to error, due to unwanted breaking change.
- `Tx2Gene`: Improved consistency of metadata return, ensuring `call` and
  `synonyms` are not defined.
- Renamed internal "acidGenomes" metadata key to "packageVersion", for
  consistency with conventions used in other Acid Genomics packages.

## AcidGenomes 0.2.8 (2021-03-02)

### Minor changes

- Relaxed validity checks for `EnsemblGenes` and `EnsemblTranscripts`.

## AcidGenomes 0.2.7 (2021-02-26)

### Minor changes

- `makeGRangesFromGFF`: Improved support and code coverage for handling of
  bcbio-nextgen `ref-transcripts.gtf` genome file.

## AcidGenomes 0.2.6 (2021-02-25)

### Minor changes

- Updated basejump dependency versions.
- `EntrezGeneInfo`: Improved column formatting.

## AcidGenomes 0.2.5 (2021-02-17)

### Minor changes

- Genome downloader functions (e.g. `downloadEnsemblGenome`) now return
  relative symlinks instead of absolute paths.
- Bug fix for `mapHumanOrthologs` internal join step. Now returns `humanGeneId`
  and `humanGeneName` columns instead of `hgncId` and `hgncName` columns, which
  technically were incorrect, since these map to Ensembl.

## AcidGenomes 0.2.4 (2021-02-13)

### New functions

- `EntrezGeneInfo`: New utility for obtaining gene annotations from NCBI.

### Major changes

- `geneSynonyms`: Reworked internal code, extending `EntrezGeneInfo`.

## AcidGenomes 0.2.3 (2021-02-10)

### Minor changes

- Reverted back to using vroom as importer for `HGNC` and `MGI2Ensembl`.

## AcidGenomes 0.2.2 (2021-02-09)

### Minor changes

- Reduced the number of reexported functions.

## AcidGenomes 0.2.1 (2021-02-08)

### Minor changes

- Now including some reexports from GenomicRanges and IRanges.

## AcidGenomes 0.2.0 (2021-02-02)

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
- GFF, TxDb, and ensembldb parser functions now default to
  `ignoreVersion = FALSE`. Previous releases of AcidGenomes and basejump had
  this set to `ignoreVersion = TRUE` by default. Note that both modes are now
  non-destructive.
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
