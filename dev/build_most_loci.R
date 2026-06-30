# Build alignments for ALL buildable loci (every locus EXCEPT B), all alignment
# types, and time it. HLA-B's cDNA build throws a pre-existing upstream error
# ("number of items to replace is not a multiple of replacement length") in BOTH
# the original and optimized code, which crashes alignmentFull(loci="all"). To
# measure the optimization at near-full scale we build the 40-locus set without
# B and compare original vs optimized for identical() output.
#
# Usage:  Rscript build_most_loci.R <label>   # label = "original" | "optimized"

.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/
suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))

label <- commandArgs(trailingOnly = TRUE)
label <- if (length(label)) label[1] else "run"
VERSION <- "3.64.0"
out <- sprintf(file.path(pkg_root, "dev/most_loci_%s.rds"), label)

# Every gazetteer locus across the three sources, minus the B-cDNA offender.
loci <- setdiff(unique(c(HLAgazetteer$nuc, HLAgazetteer$gen, HLAgazetteer$prot)), "B")

cat(sprintf("[%s] building %d loci (all except B), alignType=all, version=%s ...\n",
            label, length(loci), VERSION))
t <- system.time({
  res <- alignmentFull(loci = loci, alignType = "all", version = VERSION)
})
saveRDS(res, out)

cat(sprintf("\n[%s] DONE\n", label))
cat("elapsed (s):", round(unname(t["elapsed"]), 1),
    " user/CPU (s):", round(unname(t["user.self"]), 1), "\n")
cat("loci per type: ",
    paste(sprintf("%s=%d", names(res)[1:4],
                  vapply(res[1:4], length, integer(1))), collapse="  "), "\n")
cat("object size (MB):", round(as.numeric(object.size(res)) / 1024^2, 1), "\n")
cat("saved ->", out, "\n")
