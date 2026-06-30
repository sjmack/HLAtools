# Round-trip invariant tests across the GL String <-> UNIFORMAT conversions.
# These guard the conversion logic for cases where the transformation is known
# to be an identity round trip. NOTE: not every GL/UNIFORMAT string round-trips
# identically (cartesian-product expansion in UNtoGL can reorder/restructure
# ambiguities), so invariants are asserted only where they currently hold.

test_that("GL -> UNIFORMAT -> GL is invariant for a simple ambiguity string", {
  gl <- "HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02"
  uni <- suppressMessages(GLStoUNI(gl))
  back <- suppressMessages(UNItoGLS(uni))
  expect_identical(back, gl)
})

test_that("UNIFORMAT -> GL -> UNIFORMAT is invariant for a simple ambiguity string", {
  uni <- "A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02"
  gl <- suppressMessages(UNItoGLS(uni, pre = FALSE))
  back <- suppressMessages(GLStoUNI(gl))
  expect_identical(back, uni)
})

test_that("GLtoUN and UNtoGL are mutual inverses on the canonical example", {
  gl <- "HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02"
  uni <- "A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02"
  expect_identical(GLtoUN(gl), uni)
  expect_identical(UNtoGL(uni), gl)
})

test_that("multiUNItoGLS followed by multiGLStoUNI is structurally valid", {
  # Round-tripping the whole bundled example is NOT identity (cartesian
  # expansion restructures ambiguities), so we assert that the forward and
  # reverse passes both run and produce same-length character vectors rather
  # than asserting exact equality with the original.
  uvec <- unname(as.vector(UNIFORMAT.example[2]))[[1]][1:5]
  gvec <- suppressMessages(multiUNItoGLS(uvec))
  uvec_back <- suppressMessages(multiGLStoUNI(gvec))
  expect_type(gvec, "character")
  expect_type(uvec_back, "character")
  expect_equal(length(gvec), length(uvec))
  expect_equal(length(uvec_back), length(uvec))
})

test_that("multiGLStoUNI on the bundled GL Strings yields valid UNIFORMAT", {
  out <- suppressMessages(multiGLStoUNI(GLstring.ex[[2]][1:10]))
  # All converted strings should pass UNIFORMAT validation.
  expect_true(all(vapply(out,
                         function(x) suppressMessages(validateUniformat(x)),
                         logical(1))))
})
