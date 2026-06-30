# Characterization tests for R/GLupdater.R version/allele/GL-String translation
# helpers: redec(), GLV(), GLV2(), GLVhelper(), GLvalidate(), translateAllele(),
# translateGLstring(), updateGL(). All use the bundled alleleListHistory; only
# versions present locally (<= 3.59.0) are used so no network fallback fires.

test_that("redec reinserts version decimals from column-name form", {
  # Documented @example: redec("X3090") -> "3.09.0".
  expect_identical(redec("X3090"), "3.09.0")
  expect_identical(redec("X300"), "3.0.0")
  expect_identical(redec("X25"), "2.5")
})

test_that("GLV extracts the column-form version from a GL String Code", {
  # Documented @example: GLV("hla#3.25.0#HLA-B15:35") -> "X3250".
  expect_identical(GLV("hla#3.25.0#HLA-B15:35"), "X3250")
})

test_that("GLV2 formats a dot-delimited version into column form", {
  # Documented @examples.
  expect_identical(GLV2("3.34.0"), "X3340")
  # A 3-element version like 3.0.0 has a zero inserted -> X3000.
  expect_identical(GLV2("3.0.0"), "X3000")
})

test_that("GLVhelper returns candidate releases for truncated versions", {
  # "2.25" expands to the available 2.25.x releases (order is current behavior).
  expect_setequal(suppressMessages(GLVhelper("2.25")), c("2.25.0", "2.25.1", "2.25.2"))
  # "3.9.0" resolves to the canonical "3.09.0".
  expect_identical(suppressMessages(GLVhelper("3.9.0")), "3.09.0")
})

test_that("GLvalidate normalizes a well-formed GL String Code", {
  # Documented @example: namespace 'ha' is corrected to 'hla'; hla- -> HLA-.
  expect_identical(suppressMessages(GLvalidate("ha#3.25.0#hla-B15:35")),
                   "hla#3.25.0#HLA-B15:35")
})

test_that("GLvalidate rejects malformed GL String Codes", {
  expect_false(suppressMessages(GLvalidate("a#b#c#d")))   # too many '#'
  expect_false(suppressMessages(GLvalidate("a#b")))       # too few '#'
  expect_message(GLvalidate("a#b#c#d"), "too many")
  expect_message(GLvalidate("a#b"), "too few")
})

test_that("translateAllele translates single alleles across local versions", {
  # Documented @example: A*01:01 from 3.01.0 to 2.20.0 -> "A*0101".
  expect_identical(suppressMessages(translateAllele("A*01:01", "3.01.0", "2.20.0", FALSE)),
                   "A*0101")
  # expand = TRUE returns a slash-delimited list of all matching full names.
  exp <- suppressMessages(translateAllele("A*0101", "2.20.0", "3.01.0", TRUE))
  expect_identical(
    exp,
    "A*01:01:01:01/A*01:01:01:02N/A*01:01:02/A*01:01:03/A*01:01:04/A*01:01:05/A*01:01:06/A*01:01:07/A*01:01:08/A*01:01:09/A*01:01:10/A*01:01:11/A*01:01:12/A*01:01:13/A*01:01:14/A*01:01:15/A*01:01:16/A*01:01:17/A*01:01:18/A*01:01:19"
  )
  # Same-version expansion lists all sub-alleles.
  expect_identical(
    suppressMessages(translateAllele("B*57:01", "3.12.0", "3.12.0", TRUE)),
    "B*57:01:01/B*57:01:02/B*57:01:03/B*57:01:04/B*57:01:05/B*57:01:06/B*57:01:07/B*57:01:08/B*57:01:09/B*57:01:10/B*57:01:11/B*57:01:12/B*57:01:13/B*57:01:14"
  )
})

test_that("translateAllele errors on versions not loaded locally", {
  expect_error(translateAllele("A*01:01", "9.99.0", "2.20.0"),
               "is not loaded in the HLAtools package")
})

test_that("translateGLstring translates a full GL String across versions", {
  out <- suppressMessages(translateGLstring(GLstring.ex$Gl.String[1], "3.01.0", "2.15.0"))
  expect_identical(
    out,
    "HLA-A*0201~HLA-Cw*0702~HLA-B*0702~HLA-DRB1*1501~HLA-DQB1*0602~HLA-DPB1*0402+HLA-A*0101~HLA-Cw*0602~HLA-B*5701~HLA-DRB1*0701~HLA-DQB1*0303~HLA-DPB1*0401"
  )
  # An invalid GL String (with a forbidden character) returns FALSE.
  expect_false(suppressMessages(translateGLstring("A^B|C?D@E", "3.01.0", "2.15.0")))
})

test_that("updateGL updates a GL String Code, prepending the namespace and 'to' version", {
  out <- suppressMessages(updateGL(GLSC.ex$GL.String.Code[1], "2.15.0", FALSE, FALSE))
  expect_identical(
    out,
    "hla#2.15.0#HLA-A*0201~HLA-Cw*0702~HLA-B*0702~HLA-DRB1*1501~HLA-DQB1*0602~HLA-DPB1*0402+HLA-A*0101~HLA-Cw*0602~HLA-B*5701~HLA-DRB1*0701~HLA-DQB1*0303~HLA-DPB1*0401"
  )
  # The result carries the requested 'to' version in its second field.
  expect_true(startsWith(out, "hla#2.15.0#"))
})

test_that("updateGL errors when 'to' is not a loaded version", {
  expect_error(updateGL(GLSC.ex$GL.String.Code[1], "9.99.0"),
               "is not loaded in the HLAtools package")
})
