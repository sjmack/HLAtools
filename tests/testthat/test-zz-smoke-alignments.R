# Smoke test: confirm the alignment fixture is wired up and the consumers can
# read it offline. If this fails, the per-function tests below cannot run.

test_that("HLAalignments fixture is loaded and visible", {
  expect_true(exists("HLAalignments", envir = globalenv()))
  expect_equal(HLAalignments$version, "3.64.0")
  expect_setequal(names(HLAalignments$prot), c("A", "DRB1", "DPB1"))
})

test_that("a consumer function reads the fixture", {
  # variantTable is a small, pure consumer of HLAalignments.
  vt <- variantTable("prot", "DRB1", 86)
  expect_s3_class(vt, "data.frame")
  expect_named(vt, c("Variant", "Count", "Frequency"))
})
