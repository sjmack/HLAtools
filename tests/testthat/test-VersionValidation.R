# Characterization tests for offline parts of R/VersionValidation.R:
# squashVersion(), expandVersion(), repoVersion(), checkVersion().
# validateVersion() and getLatestVersion() are NOT tested here because they reach
# out to the network when a version is not present locally.

test_that("squashVersion removes dot delimiters", {
  # Documented @example: squashVersion("3.45.0", TRUE) -> 3450 (numeric).
  expect_equal(squashVersion("3.45.0", TRUE), 3450)
  expect_identical(squashVersion("3.45.0"), "3450")
  # A numeric input is coerced to character form (no dots to strip).
  expect_identical(squashVersion(3450), "3450")
})

test_that("expandVersion reinserts dot delimiters", {
  # Documented @example: expandVersion(3450) -> "3.45.0".
  expect_identical(expandVersion(3450), "3.45.0")
  expect_identical(expandVersion("3450"), "3.45.0")
})

test_that("repoVersion converts to GitHub-branch form", {
  # Documented @example: repoVersion("3.05.0") -> "350" (3.00.0-3.09.0 collapse).
  expect_identical(repoVersion("3.05.0"), "350")
  # A non-zero second field is left as the 4-digit squashed form.
  expect_identical(repoVersion("3.25.0"), "3250")
})

test_that("checkVersion detects locally bundled release versions", {
  # 3.25.0 is present in the bundled alleleListHistory.
  expect_true(checkVersion("3.25.0"))
  # An absurd version absent from the bundle returns FALSE (no network call,
  # because checkVersion only inspects the local object).
  expect_false(checkVersion("9.99.0"))
})
