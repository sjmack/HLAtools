# HLAtools Performance Optimization — `alignmentFull()` / `buildAlignments()`

This document records a performance-optimization pass over `HLAtools`, focused on
`alignmentFull()` and the function it drives, `buildAlignments()`. **Every
performance change preserves behavior exactly**: the optimized code reproduces the
original output byte-for-byte (verified with `identical()` against a saved
baseline). A separate, clearly-scoped bug fix ([§10](#10-bug-fix-hla-b-cdna-crash-in-the-all-loci-build))
makes the default `alignmentFull(loci = "all")` work again.

- **Upstream:** `sjmack/HLAtools` (v1.8.1, commit `edc9cdc`)
- **Fork / working copy:** `k96nb01/HLAtools`, branch `perf-optimization`
- **Headline result:** the representative 3-locus, all-types build dropped from
  **~139 s to ~48 s (≈2.9× faster)**; at full scale (40 loci) **401 s → 159 s
  (2.5× wall, 2.9× CPU)**, output `identical()`. The performance work changes no
  output.
- **Bonus fix:** `alignmentFull(loci = "all")` previously crashed on a
  pre-existing HLA-B cDNA bug; it is fixed here and now builds the full gazetteer.

---

## 1. Motivation

`alignmentFull()` builds the `HLAalignments` object that nearly every other
alignment-aware function in the package depends on. It is slow enough that
building all loci takes several minutes, which is a recurring friction point in
day-to-day use. The goal was to make it substantially faster **without changing
any output**, so existing analyses remain valid.

---

## 2. Method: correctness first, measure before cutting

The package shipped with **no `tests/` directory at all**, so the first job was
to build a safety net, not to change code.

1. **Characterization baseline.** Ran the *original* `alignmentFull()` on a
   representative subset of loci and saved the exact result. Loci chosen to
   exercise the distinct code paths:
   - `A` — class I (prot/nuc/gen), and the redundant-cDNA download path
   - `DRB1` — the special DR-locus handling, and a very wide gDNA alignment
     (662 × 18 527)
   - `DPB1` — DP naming conventions

   The version is **pinned to a concrete release (3.64.0)**, not `"Latest"`, so
   the baseline is reproducible (IPD-IMGT/HLA release branches are immutable
   upstream).

2. **`identical()` gate after every edit.** Two harnesses (see
   [§7](#7-reproducing-the-work)) rebuild the loci and assert the new output is
   `identical()` to the baseline. A single non-identical byte fails the gate.
   This was confirmed to report *IDENTICAL* on the unmodified code first, which
   also proves `alignmentFull()` is deterministic (a prerequisite for this
   approach).

3. **Profile, don't guess.** Line-level profiling (`Rprof(line.profiling=TRUE)`)
   located the real hotspots. This step changed the strategy entirely — see
   [§4](#4-the-key-diagnosis-cpu-bound-not-network-bound). The harnesses are in
   `dev/` — see [§8](#8-reproducing-the-work).

4. **A regression test suite.** Alongside the refactor, a full `testthat`
   (edition 3) suite was written: **18 files, 102 tests, 255 expectations, all
   passing** — offline golden tests for the pure functions plus
   alignment-consumer tests run against the saved fixture
   ([§7](#7-test-suite)).

---

## 3. Environment

- R 4.6.0 on Windows 11.
- Added dev/test dependencies: `testthat`, `pkgload`, `profvis`, `bench`, plus
  the two package `Imports` that were not yet installed (`DescTools`, `fmsb`).
- Network access to `raw.githubusercontent.com/ANHIG/IMGTHLA` confirmed
  (required for any real build).

---

## 4. The key diagnosis: CPU-bound, not network-bound

A static read of the call graph suggested the cost was network I/O — a full run
issues ~135 sequential HTTP downloads. **Profiling proved otherwise.** The
baseline build reported:

```
   user  system elapsed
 124.36    7.03  139.11
```

`user` (CPU) is **124 s of 139 s elapsed** — roughly **90% of the time is CPU**,
spent parsing and reshaping the downloaded alignment text. Network was only
~10%. This redirected all effort to the in-memory parsing in `buildAlignments()`.

Line-level profiling of three representative builds (`A/cDNA`, `A/gDNA`,
`DRB1/gDNA`) found two dominant lines:

| Line | Share of CPU | What it did |
|------|-------------:|-------------|
| `buildAlignments.R` reference-distribution loop | **~35%** | per-column `for` loop over up to ~18.5k columns, re-indexing the data frame each pass |
| `buildAlignments.R` sequence split | **~25%** | `sapply(seqs, strsplit, split="*")` splitting every sequence into characters |

The download lines were ~4% combined.

---

## 5. The optimizations

All changes are commented inline in the source explaining the rationale. Each was
verified `identical()` to the baseline before moving on.

### 5.1 Vectorize the reference-distribution loops *(biggest single win)*

**File:** `R/buildAlignments.R` (the two "distributes reference sequence from
row 1" blocks).

The original distributed the reference sequence into every `"-"` cell one column
at a time:

```r
for(x in 5:ncol(DNAalignments[[loci[i]]])) {
  DNAalignments[[loci[i]]][,x][which(DNAalignments[[loci[i]]][,x]=="-")] <-
    DNAalignments[[loci[i]]][,x][1]}
```

For a gDNA alignment with ~18.5k position columns this loops 18.5k times, and
each pass re-resolves the nested `[[loci[i]]]` list lookup and a data-frame
column extraction. Replaced with a single vectorized matrix operation:

```r
dcols <- 5:ncol(DNAalignments[[loci[i]]])
dm <- as.matrix(DNAalignments[[loci[i]]][, dcols, drop = FALSE])
ddash <- which(dm == "-", arr.ind = TRUE)            # every dash, as (row, col)
if(nrow(ddash) > 0) { dm[ddash] <- dm[1, ddash[, "col"]] }  # row-1 value per column
DNAalignments[[loci[i]]][, dcols] <- dm
```

This was **~35% of CPU**. The same fix was applied to the AA/codon block.

### 5.2 `strsplit(split = "")` instead of `split = "*"` *(14× on that line)*

**File:** `R/buildAlignments.R` (the per-character sequence split).

The code split sequences into characters with `strsplit(x, split="*")`. `"*"` is
a **degenerate regex** — a quantifier with nothing to repeat — which forces
`strsplit()` onto its slow regex engine. The empty pattern `""` produces the
*identical* character split but takes `strsplit()`'s special-cased fast path.
Benchmarked on a realistic sequence vector:

| split | median | identical to `"*"`? |
|-------|-------:|:-------------------:|
| `"*"` (original) | 1.37 s | — |
| `""` (new) | **0.098 s** | **yes** |

A ~**14×** speedup on a line that was ~25–31% of CPU, by changing one character.
(The surrounding `sapply` was also collapsed to a single vectorized `strsplit`
call, since `strsplit` already vectorizes over its input.)

### 5.3 `match()`-vectorize the INDEL / EXONB labeling

**File:** `R/buildAlignments.R` (the inDel and exon-boundary inclusion blocks).

The original scanned the entire correspondence-table label row once per indel —
`O(nIndels × ncol)`:

```r
for(o in 1:length(inDels[[loci[i]]])){
  corr_table[[loci[i]]][2,][inDels[[loci[i]]][[o]]==corr_table[[loci[i]]][1,]] <- paste("INDEL", o, sep="-")
  ...
}
```

The position labels in `corr_table[1,]` were **verified to be unique** (across
cDNA and gDNA builds, including the 18 523-column DRB1 gDNA case), so `match()`
locates every indel column in a single pass with identical results:

```r
indelIdx <- match(inDels[[loci[i]]], corr_table[[loci[i]]][1,])
indelLab <- paste("INDEL", seq_along(indelIdx), sep="-")
corr_table[[loci[i]]][2, indelIdx] <- indelLab
corr_table[[loci[i]]][3, indelIdx] <- indelLab
```

### 5.4 `which(is.na())`-vectorize the position-fill loops

**File:** `R/buildAlignments.R` (the two "pastes alignment_positions into
corr_table" loops).

Scalar loops walked every column, assigning the next position value to each
non-indel cell with a manual counter. Since the `NA` cells, in column order,
simply receive `alignment_positions[1..k]`, this is one vectorized assignment:

```r
naIdx2 <- which(is.na(corr_table[[loci[i]]][2,]))
corr_table[[loci[i]]][2, naIdx2] <- alignment_positions[[loci[i]]][seq_along(naIdx2)]
```

### 5.5 Eliminate the duplicate cDNA build *(structural)*

**File:** `R/alignmentFull.R`.

`buildAlignments(locus, "cDNA")` returns **both** the codon and cDNA tables in a
single call. But `alignmentFull()` called it **twice per locus** — once for the
`nuc` list (keeping the cDNA table) and once for the `codon` list (keeping the
codon table) — re-downloading and re-parsing the identical `_nuc.txt` file and
discarding half of each result. (This was visible in the baseline: the `codon`
and `nuc` tables have identical dimensions because they come from the same
file.) The fix builds each required cDNA alignment **once** into a small cache
and reads both tables from it:

```r
cdnaLoci <- unique(c(as.character(NL1), as.character(NL4)))
cdnaLoci <- cdnaLoci[!is.na(cdnaLoci) & nzchar(cdnaLoci)]
cdnaCache <- vector(mode = 'list', length = length(cdnaLoci))
names(cdnaCache) <- cdnaLoci
for(lc in cdnaLoci){
  cdnaCache[[lc]] <- suppressWarnings(buildAlignments(lc, "cDNA", version = version)[[1]])
}
# cList  (nuc)   reads cdnaCache[[locus]][2]
# codonList (codon) reads cdnaCache[[locus]][1]
```

On a full run this removes roughly a quarter of *all* `buildAlignments` calls —
both a download and a full re-parse per locus.

---

## 6. Benchmark results

All builds pinned to IPD-IMGT/HLA release 3.64.0.

**Full `alignmentFull(loci = c("A","DRB1","DPB1"), alignType = "all")`:**

| | elapsed | vs original |
|--|--------:|------------:|
| Original (baseline) | 139.1 s | 1.0× |
| Optimized | ~48 s | **≈2.9×** |

(Wall-clock now includes a larger relative share of network time, so it varies
run to run; CPU is the stabler metric below.)

**CPU work — profiler sampled time, three representative `buildAlignments`
calls (`A/cDNA`, `A/gDNA`, `DRB1/gDNA`):**

| Stage | sampled CPU |
|-------|------------:|
| Original | 50.75 s |
| After §5.1 + §5.2 | 42.31 s |
| Final (all changes) | **18.47 s** |

≈**2.75× less CPU**, and the profile is now *flat* — no single line dominates.

**Per-`buildAlignments` harness (5 calls covering every hotspot):** 88.6 s → 37.9 s
(**2.34×**).

### Full-scale: 40-locus build

To confirm the result at near-full scale, the entire gazetteer was built — every
locus across all alignment types — **except HLA-B**, whose cDNA build crashes in
the *original* code (a pre-existing upstream bug now fixed in this fork; see
[§10](#10-bug-fix-hla-b-cdna-crash-in-the-all-loci-build)). The original cannot
build B at all, so the apples-to-apples comparison uses the 40-locus set
(prot = 26, codon = 40, nuc = 40, gen = 40; ~1.96 GB object):

| | elapsed | CPU (user) |
|--|--------:|-----------:|
| Original | 401.3 s | 336.2 s |
| Optimized | **159.4 s** | **116.5 s** |
| **Speedup** | **2.52×** | **2.89×** |

Critically, the two 40-locus results were compared with `identical()` and are
**identical across every locus and alignment type** — full-scale proof that the
optimization changes nothing about the output.

---

## 7. Test suite

The package had no tests. The new `testthat` suite (`tests/testthat/`) is the
regression net:

- **Offline golden tests** (12 files, ~160 expectations) for the pure / data-only
  functions: `alleleTrim`, `getField`, GL ↔ UNIFORMAT conversion and
  round-trips, version utilities, locus/motif validators, `BDtoPyPop`, and more.
- **Alignment-consumer tests** (6 files, ~95 expectations) for the functions that
  read `HLAalignments` — `compareSequences`, `motifMatch`, the
  `alignmentSearch`/`customAlign` family, `posSort`, `validateAllele` — run
  against the saved baseline fixture so they need no network.

**Total: 102 tests / 255 expectations, all passing.**

These tests also documented several pre-existing behavior quirks, captured as
characterization (current behavior is locked, not "fixed"):

- `BDtoPyPop()` names its second returned list element `...neagtive` (misspelled).
- A `multiGLStoUNI()` roxygen `@example` passes a non-existent `version` object.
- `posSort()` returns its result in the input's type (numeric in → numeric out),
  not always character as the docs imply.
- `validateAllele()` on an invalid locus *errors* (inside `checkSource`) rather
  than returning `FALSE`.

---

## 8. Reproducing the work

The `dev/` folder (build-ignored) holds the harnesses:

| Script | Purpose |
|--------|---------|
| `dev/build_baseline.R` | Build & save the `alignmentFull` characterization baseline + fixture |
| `dev/verify_baseline.R` | Rebuild the baseline loci and assert `identical()` (the milestone gate, ~48 s) |
| `dev/verify_fast.R` | Faster per-`buildAlignments` `identical()` check for rapid iteration |
| `dev/profile.R` | Line-level `Rprof` profile of `buildAlignments` |
| `dev/run_tests.R` | Run the full `testthat` suite via `pkgload::load_all` |
| `dev/build_most_loci.R` | Build the 40-locus (all-except-B) set; `<label>` selects the output file |
| `dev/compare_most_loci.R` | Assert the original vs optimized 40-locus builds are `identical()` |
| `dev/diagnose_loci.R` | Per-locus build check that isolated the HLA-B cDNA failure |
| `dev/validate_B_fix.R` | Confirm HLA-B builds for all sources after the fix |
| `dev/verify_fix_noop.R` | Prove the B fix leaves all 40 other loci `identical()` |
| `dev/imgt_release_check.R` | Smoke-test all loci against an IMGT release (used by CI; see [§11](#11-continuous-imgt-release-monitoring)) |

Each `dev/` script resolves the repository root from its own file location, so it
has no hardcoded paths and runs from any working directory or in CI.

The before/after at full scale was produced by checking out the original two R
files (`git checkout <upstream> -- R/alignmentFull.R R/buildAlignments.R`),
running `build_most_loci.R original`, restoring the optimized files, running
`build_most_loci.R optimized`, then `compare_most_loci.R`.

Run any of them with, e.g. (each script locates the repo root from its own
path, so it works from any working directory):

```pwsh
Rscript dev/verify_baseline.R
```

Fixtures live in `tests/testthat/fixtures/` (`alignmentFull_baseline.rds`; the
per-call fast baseline is regenerated on demand).

---

## 9. What was *not* changed (and why)

After the changes above, the profile is flat. The remaining costs are spread
thin and were left alone as poor risk/reward:

- **Final list-assembly copies** (`c(final_alignment[[...]], <data frame>)`) —
  inherent to building the nested result structure; hard to change without
  risking output differences.
- **Whitespace normalization / `str_squish`** on the raw downloaded lines —
  string parsing where any change risks altering output.
- **Parallel / cached downloads** — network is now only ~10% of the time, so the
  payoff is small relative to the complexity (and concurrent requests to the
  ANHIG repo are best avoided). This is the natural next lever if a full-repo
  build is still too slow.

---

## 10. Bug fix: HLA-B cDNA crash in the all-loci build

While running the full-scale test, the **default `alignmentFull(loci = "all")`
was found to crash** on release 3.64.0 with:

```
Error in corr_table[[loci[i]]][1, ] <- names(HLAalignments[[loci[i]]][5:ncol(...)]) :
  number of items to replace is not a multiple of replacement length
```

This is **a pre-existing bug in the upstream code**, *not* introduced by this
optimization work — the failing line is original code that was not modified, and
the error reproduces identically on the unmodified `sjmack/HLAtools` v1.8.1.

**This fork fixes it.**

**Root cause.** A per-locus diagnostic (`dev/diagnose_loci.R`) isolated the
failure to **exactly one build: HLA-B `cDNA`** (HLA-B's `gDNA` and `AA` build
fine, as do all other loci). In `B_nuc.txt`, a handful of alleles — e.g.
`B*44:568Q` and `B*51:197` — appear in a **trailing alignment block that the
reference allele `B*07:02:01:01` does not span** (the reference spans 17 blocks;
these alleles span 18). This is a genuine IMGT data structure: those alleles
carry 3′ sequence extending past the reference, which leaves their assembled
sequences a different length from the reference (1738 and 1458 vs the reference's
1473). `do.call(rbind, ...)` — with its warning suppressed — then silently
recycled these ragged rows to the longest one, corrupting the alignment width and
throwing the cryptic error at the correspondence-table fill (whose width is the
reference length).

**The fix** (`R/buildAlignments.R`, just after the per-character split). The
reference allele defines the position coordinate system, so every cDNA/gDNA
sequence is normalized to the reference length before the bind: sequence
extending **beyond** the reference is trimmed (it has no reference-relative
position), and **short** partial sequences are padded with `"."` — the same fill
the function already uses for amino-acid premature termination. The reference row
is never altered, so the indel/exon-boundary detection that reads row 1 is
unaffected. A `message()` reports the locus, source, and number of sequences
normalized.

**Validation.**
- `alignmentFull(loci = "all")` now **completes** for the full gazetteer
  (prot = 27, codon = 41, nuc = 41, gen = 41).
- Across that full build, **only the two HLA-B cDNA alleles triggered the
  normalization** — every other locus is a no-op.
- The 40 non-B loci were compared table-by-table (**146 tables**) against the
  pre-fix build and are **identical**; the A/DRB1/DPB1 baseline is also unchanged.
  The fix is surgical — it makes HLA-B buildable and changes nothing else.
- A network-guarded regression test
  (`tests/testthat/test-buildAlignments-B-regression.R`) builds HLA-B cDNA and
  asserts the table is rectangular (11 110 × 1 477) with both formerly-ragged
  alleles present.

**Note for upstream.** This is a behavior change for the two affected HLA-B cDNA
alleles only — they are normalized to the reference frame. Because the original
code could not produce *any* HLA-B cDNA output (it crashed), there is no prior
baseline for those alleles to preserve. The full-scale benchmark in
[§6](#6-benchmark-results) still reports the 40-locus (ex-B) set so the before/
after is an apples-to-apples comparison against the original, which cannot build
B at all.

---

## 11. Continuous IMGT release monitoring

The package downloads alignment data live from the ANHIG/IMGTHLA GitHub
repository, which IMGT updates roughly quarterly. As the HLA-B case showed, a new
release can introduce data that the alignment parser chokes on. To catch that
proactively, this fork adds a scheduled GitHub Actions workflow
(`.github/workflows/imgt-release-check.yml`) that builds the alignments against
the live release data and turns red if anything breaks.

The workflow runs `dev/imgt_release_check.R`, which:

1. Runs the real user path — `alignmentFull(loci = "all", alignType = "all")` —
   against `"Latest"`.
2. **On success:** reports the version and locus counts, exits 0.
3. **On failure:** runs a per-locus / per-source sweep and prints exactly which
   build(s) broke (the same diagnostic that isolated the HLA-B bug), then exits
   non-zero so the job fails and the maintainer is notified.

It also runs the unit + regression test suite (with `NOT_CRAN` set so the
network-guarded HLA-B regression test runs against the live data).

**Triggers.** Weekly (`cron`), plus a manual *Run workflow* button that accepts a
specific version to test (e.g. the day a new release drops). Note that GitHub
only fires *scheduled* workflows from a repository's **default branch**, so the
weekly schedule begins once this work is on the default branch; the manual
trigger works from any branch immediately.

**Scope.** The loci tested come from the bundled `HLAgazetteer`, so this verifies
that *existing* loci still build under a new release. A release that adds an
entirely new locus would also require rebuilding the gazetteer
(`buildGazetteer()` / `updateAll()`) — a separate concern from this check.

---

*Generated as part of the HLAtools performance-optimization effort. All
performance optimizations are verified to produce output `identical()` to
HLAtools v1.8.1; the HLA-B fix is the one intentional, validated behavior change.*
