# Lessons

## "NCBI Gene ID" and "Entrez ID" are the same thing

"NCBI Gene ID" and "Entrez ID" refer to the same identifiers. NCBI Gene is the
database; Entrez is the overarching retrieval system that once named the database
"Entrez Gene." Both terms are still valid in bioinformatics. The codebase uses
`ncbiGeneId` as the canonical column name (not `entrezId`) to match the official
database name, but either term is acceptable in conversation and comments.

## Verification standard: pass `AcidDevTools::check()`

Never mark R package work complete without running `AcidDevTools::check()`.

`testthat::test_local()` alone is insufficient — it skips R CMD check, lintr,
documentation validation, and other checks that `AcidDevTools::check()` runs.

Always run `AcidDevTools::check()` as the final verification step before
declaring any code change done.
