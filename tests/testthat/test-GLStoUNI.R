# Characterization tests for R/GLStoUNI.R: validateGLstring(), GLtoUN(),
# GLStoUNI(), multiGLStoUNI(). Uses the bundled GLstring.ex dataset.

test_that("validateGLstring enforces version-specific character sets", {
  # Documented @example: a well-formed v1.0 GL String validates.
  expect_true(validateGLstring("HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02", "1.0"))
  # The "?" operator is forbidden in v1.0 but allowed in v1.1.
  expect_false(suppressMessages(validateGLstring("A?B", "1.0")))
  expect_true(validateGLstring("A?B", "1.1"))
  # An unsupported version string returns FALSE.
  expect_false(validateGLstring("A+B", "2.0"))
  # Forbidden characters emit a message.
  expect_message(validateGLstring("A?B", "1.0"), "not permitted in GL String")
})

test_that("GLtoUN translates GL Strings to UNIFORMAT", {
  # Documented @example 1.
  expect_identical(
    GLtoUN("HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02"),
    "A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02"
  )
  # Documented @example 2 (mixed operators).
  expect_identical(
    GLtoUN("A+B/C~D^G|E^W+X/Y^Z+J"),
    "A,B|A,C~D G|E W,X|W,Y Z,J"
  )
})

test_that("GLStoUNI wraps GLtoUN and rejects invalid strings", {
  expect_identical(
    suppressMessages(GLStoUNI("HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02")),
    "A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02"
  )
  expect_identical(
    suppressMessages(GLStoUNI("A+B/C~D^G|E^W+X/Y^Z+J")),
    "A,B|A,C~D G|E W,X|W,Y Z,J"
  )
  # The "?" operator has no UNIFORMAT cognate; GLStoUNI returns FALSE.
  expect_false(suppressMessages(GLStoUNI("A?B")))
})

test_that("multiGLStoUNI converts a vector of GL Strings", {
  # NOTE: the roxygen @example passes a non-existent 'version' as the 2nd
  # positional arg (which is actually 'prefix'); calling without it is the
  # intended usage. The first 5 GL Strings convert to UNIFORMAT as below.
  out <- suppressMessages(multiGLStoUNI(GLstring.ex[[2]][1:5]))
  expect_type(out, "character")
  expect_equal(length(out), 5L)
  expect_identical(
    out[1],
    "A*02:01~C*07:02~B*07:02~DRB1*15:01~DQB1*06:02~DPB1*04:02,A*01:01~C*06:02~B*57:01~DRB1*07:01~DQB1*03:03~DPB1*04:01"
  )
})

test_that("multiGLStoUNI converts a data frame and preserves its identifier column", {
  out <- suppressMessages(multiGLStoUNI(GLstring.ex[1:5, 1:2]))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 5L)
  # The identifier column (first column) is untouched.
  expect_equal(out$Relation, rep("Subject", 5L))
  # The GL String column has been translated to UNIFORMAT (commas present).
  expect_true(all(grepl(",", out$Gl.String)))
})
