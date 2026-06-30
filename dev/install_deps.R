# Install packages needed for the optimization project:
#  - DescTools, fmsb : package Imports not yet installed (needed for load_all)
#  - pkgload         : load_all() to load the package from source
#  - testthat        : the test suite
#  - profvis, bench  : profiling and benchmarking the optimization work
need <- c("DescTools","fmsb","pkgload","testthat","profvis","bench")
inst <- rownames(installed.packages())
need <- need[!need %in% inst]
if (length(need)) {
  install.packages(need, repos = "https://cloud.r-project.org", quiet = FALSE)
} else cat("nothing to install\n")
# Report final status
inst2 <- rownames(installed.packages())
cat("---FINAL---\n")
for (p in c("DescTools","fmsb","pkgload","testthat","profvis","bench"))
  cat(sprintf("%-12s %s\n", p, if (p %in% inst2) "OK" else "STILL MISSING"))
