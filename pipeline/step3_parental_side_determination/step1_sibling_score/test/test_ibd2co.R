## Unit tests for ibd2co (one-group THORIN labels, v1.2.0 release and v1.3).
## Run from pipeline/step3_parental_side_determination/step1_sibling_score/:
##   Rscript test/test_ibd2co.R
## Exits with an error if any test fails.

source('src/utils.R')

BP_PER_CM <- 1e6
TOL_BP    <- 1e4   # crossover position must be within 10 kb of the truth

## Build a THORIN-like .ibd table from two true sharing tracks.
## h0 / h1: data.frame(from_cM, to_cM, shared) covering the same chromosome.
## Labels: A = only h0 shared, B = only h1, C = neither, D = both.
make_ibd <- function(h0, h1, target='T1', chr='chr1'){
  cuts <- sort(unique(c(h0$from_cM, h0$to_cM, h1$from_cM, h1$to_cM)))
  out <- data.frame()
  for (i in seq_len(length(cuts) - 1)){
    a <- cuts[i]; b <- cuts[i + 1]; mid <- (a + b) / 2
    s0 <- h0$shared[h0$from_cM <= mid & h0$to_cM > mid]
    s1 <- h1$shared[h1$from_cM <= mid & h1$to_cM > mid]
    lab <- if (s0 & !s1) 'A' else if (!s0 & s1) 'B' else if (!s0 & !s1) 'C' else 'D'
    out <- rbind(out, data.frame(CHR=chr, start=a * BP_PER_CM, end=b * BP_PER_CM - 1000,
                                 Prob=lab, length_CM=b - a, target=target, UPD='U'))
  }
  ## THORIN writes one row per run of identical labels: merge identical neighbours
  r <- rle(out$Prob); e <- cumsum(r$lengths); s <- c(1, head(e, -1) + 1)
  data.frame(CHR=chr, start=out$start[s], end=out$end[e], Prob=r$values,
             length_CM=mapply(function(x, y) sum(out$length_CM[x:y]), s, e),
             target=target, UPD='U')
}

check <- function(res, expected, name, verbose=TRUE){
  ## expected: data.frame(cM, haplotype)
  ok <- !is.null(res) && nrow(res) == nrow(expected)
  if (ok){
    res <- res[order(res$haplotype, res$bp), ]
    expected <- expected[order(expected$haplotype, expected$cM), ]
    ok <- all(res$haplotype == expected$haplotype) &&
          all(abs(res$bp - expected$cM * BP_PER_CM) <= TOL_BP)
  }
  if (!ok){
    cat('FAILED:', name, '\nGot:\n'); print(res); cat('Expected:\n'); print(expected)
    stop('test failed: ', name)
  }
  if (verbose) cat('ok  -', name, '\n')
}

## Test 1: A -> D -> B -> C, with a short (1 cM) C flicker inside the A run.
## True crossovers: h1 at 30, h0 at 50, h1 at 80.
h0 <- data.frame(from_cM=c(0, 10, 11, 50), to_cM=c(10, 11, 50, 100), shared=c(TRUE, FALSE, TRUE, FALSE))
h1 <- data.frame(from_cM=c(0, 30, 80),     to_cM=c(30, 80, 100),    shared=c(FALSE, TRUE, FALSE))
ibd <- make_ibd(h0, h1)
stopifnot(identical(ibd$Prob, c('A', 'C', 'A', 'D', 'B', 'C')))
check(ibd2co(ibd, 'T1', 3, 1),
      data.frame(cM=c(30, 50, 80), haplotype=c('H1', 'H0', 'H1')),
      'A->D->B->C with short flicker ignored')

## The legacy function gets this case wrong (documents the bug).
old <- ibd2co_v12(ibd, 'T1', 3, 1)
stopifnot(any(old$haplotype == 'H0' & abs(old$bp - 30 * BP_PER_CM) < 2e6))
cat('ok  - legacy ibd2co_v12 puts the A->D boundary (h1 crossover) on H0, as expected\n')

## Test 2: D -> C, where both haplotypes change at the same position.
h0 <- data.frame(from_cM=c(0, 60),     to_cM=c(60, 100),     shared=c(TRUE, FALSE))
h1 <- data.frame(from_cM=c(0, 40, 60), to_cM=c(40, 60, 100), shared=c(FALSE, TRUE, FALSE))
ibd <- make_ibd(h0, h1)
stopifnot(identical(ibd$Prob, c('A', 'D', 'C')))
check(ibd2co(ibd, 'T1', 3, 1),
      data.frame(cM=c(40, 60, 60), haplotype=c('H1', 'H0', 'H1')),
      'A->D->C, double change at D->C')

## Test 3: no crossover at all, and wrong target / chromosome.
h0 <- data.frame(from_cM=0, to_cM=100, shared=TRUE)
h1 <- data.frame(from_cM=0, to_cM=100, shared=FALSE)
ibd <- make_ibd(h0, h1)
stopifnot(is.null(ibd2co(ibd, 'T1', 3, 1)))
stopifnot(is.null(ibd2co(ibd, 'OTHER', 3, 1)))
stopifnot(is.null(ibd2co(ibd, 'T1', 3, 2)))
cat('ok  - no crossover / other target / other chromosome return NULL\n')

## Test 4: random tracks. Crossovers at least 4 cM apart on each haplotype,
## so every run is >= 3 cM and all crossovers must be recovered exactly.
set.seed(1)
random_track <- function(n_co){
  repeat {
    co <- sort(runif(n_co, 4, 196))
    if (n_co == 0 || all(diff(c(0, co, 200)) >= 4)) break
  }
  b <- c(0, co, 200)
  first <- sample(c(TRUE, FALSE), 1)
  data.frame(from_cM=head(b, -1), to_cM=b[-1],
             shared=rep(c(first, !first), length.out=length(b) - 1))
}
for (i in 1:300){
  h0 <- random_track(rpois(1, 2)); h1 <- random_track(rpois(1, 2))
  ibd <- make_ibd(h0, h1, chr='chr7')
  expected <- rbind(data.frame(cM=h0$from_cM[-1], haplotype=rep('H0', nrow(h0) - 1)),
                    data.frame(cM=h1$from_cM[-1], haplotype=rep('H1', nrow(h1) - 1)))
  res <- ibd2co(ibd, 'T1', 3, 7)
  if (nrow(expected) == 0){ stopifnot(is.null(res)); next }
  check(res, expected, paste('random track', i), verbose=FALSE)
}
cat('ok  - 300 random tracks recovered exactly\n')

cat('\nAll ibd2co tests passed.\n')
