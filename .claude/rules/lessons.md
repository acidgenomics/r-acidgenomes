# Lessons

## Verification standard: pass `AcidDevTools::check()`

Never mark R package work complete without running `AcidDevTools::check()`.

`testthat::test_local()` alone is insufficient — it skips R CMD check, lintr,
documentation validation, and other checks that `AcidDevTools::check()` runs.

Always run `AcidDevTools::check()` as the final verification step before
declaring any code change done.
