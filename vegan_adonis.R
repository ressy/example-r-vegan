# Example demonstrating the vegan package's adonis() function, for mock distance
# matrices.  We should see a significant p-value for sample groups with large
# differences and the opposite for sample groups that look very similar.


# Setup -------------------------------------------------------------------


set.seed(1)

# Create a mock distance matrix with a number of separate groups of similar
# samples, with the groups separated by a given distance minimum.  I think these
# are more like dissimilarity scores since I'm doing nothing to handle the
# triangle rule here (so sample s1 and s2 could be close to s3 but far from each
# other).
mock_dist_mat <- function(samples_per_group=10, num_groups=2, min_dist=10, sd=1, lows=1) {
  N <- samples_per_group
  dists <- do.call(cbind, lapply(1:num_groups, function(i) {
    left <- if (i > 1)
      (num_groups+2-i):num_groups
    else
      c()
    right <- 1:(num_groups+1-i)
    idx <- c(left, right)
    mats <- lapply(1:num_groups, function(i) {
      m <- matrix(rnorm(N^2, sd=sd)^2, nrow=N)
      if (! i %in% lows)
        m <- m + min_dist
      m
    })
    do.call(rbind, mats[idx])
  }))
  s <- paste0("s", 1:(N*num_groups))
  rownames(dists) <- s
  colnames(dists) <- s
  # Being lazy here and just letting as.dist() take the lower diagonal to make
  # it symmetric.  The other half of the matrix is ignored.
  as.dist(dists)
}

dist_heatmap <- function(dists, ...) {
  pheatmap::pheatmap(dists,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     labels_row = attr(dists, "Labels"),
                     labels_col = attr(dists, "Labels"),
                     ...)
}


# Give a list of vectors with every permutation of the items in the given
# vector.  (How is this not in base?)
permutation <- function(rem, vec=c(), maxdepth=0, depth=0) {
  halt <- (maxdepth > 0 && depth > maxdepth)
  x <- mapply(function(r, i) {
    if (length(rem) == 1 || halt) {
      c(vec, rem)
    } else {
      permutation(rem[-i], c(vec, r), maxdepth, depth+1)
    }
  }, rem, 1:length(rem), SIMPLIFY = FALSE)
  if (length(rem) > 1 && ! halt) {
    x <- unlist(x, recursive = FALSE)
  }
  x
}


# Significant -------------------------------------------------------------


N <- 5
num_groups <- 3
min_dist <- 0.1
sd <- 0.2
dists <- mock_dist_mat(N, num_groups = num_groups, min_dist, sd)

# Sample attributes: to start with, just grouping the identical samples together.
# Note that adonis does NOT care about the row names, just order.
attrs <- data.frame(SampleGroup = rep(LETTERS[1:num_groups], each=N))
rownames(attrs) <- attr(dists, "Labels")
# So, don't do this:
#attrs <- attrs[shuffle(attrs), , drop=F]

# This creates an adonis object, with a list of objects including an anova
# object (as from anova()) with the main results.
result <- adonis(dists ~ SampleGroup, data = attrs)
result$aov.tab
# num permutations for "entire set" is factorial(attr(dists, "Size")) - 1
p <- result$aov.tab$`Pr(>F)`[1]
dist_heatmap(dists, main = paste0("Significant, p=", p))


# Insignificant -----------------------------------------------------------


min_dist <- 0
dists <- mock_dist_mat(N, num_groups = num_groups, min_dist, sd)
result <- adonis(dists ~ SampleGroup, data = attrs)
result$aov.tab
p <- result$aov.tab$`Pr(>F)`[1]
dist_heatmap(dists, main = paste0("Insignificant, p=", p))


# One Different -----------------------------------------------------------


# Here one group is different but the other two are not.
min_dist <- 0.1
dists <- mock_dist_mat(N, num_groups = num_groups, min_dist, sd, lows=1:2)
result <- adonis(dists ~ SampleGroup, data = attrs)
result$aov.tab
p <- result$aov.tab$`Pr(>F)`[1]
dist_heatmap(dists, main = paste0("One Different, p=", p))


# Stats Check -------------------------------------------------------------


# Let's see if I can get the same answers myself, just mucking through the math
# in base R.
# Referring to:
# https://en.wikipedia.org/wiki/Permutational_analysis_of_variance
# http://img2.timg.co.il/forums/1_124959686.pdf
# Here I'll re-run adonis() with explicitly-selected permutations.

# From here on I haven't been able to make things exactly match up when trying
# to define the permutations myself.  Possibly I'm misunderstanding how the
# adonis documentation defines the matrix.  (Does "each row gives the permuted
# indices" not mean what I think it does?)

nn <- N*num_groups

perms <- do.call(rbind, permutation(1:nn, maxdepth = nn))
perms <- unique(perms)

result <- adonis(dists ~ SampleGroup, data = attrs, permutations = perms)

# First off, sums of squares for SampleGroup.
# Convenience function for these summations across dist matrices.
diaggr <- function(d, func) {
  1/(nrow(d))*sum(sapply(1:(nrow(d)-1), function(i) {
    sum(sapply((i+1):nrow(d), function(j) {
        func(i, j)
      }))
  }))
}
d <- as.matrix(dists)
# total sum of squares:
sst <- diaggr(d, function(i, j) d[i, j]^2 )
# within-group sum of squares:
eps <- function(a, i, j) (a[i] == a[j]) + 0
ssw <- num_groups*diaggr(d, function(i, j) d[i, j]^2 * eps(attrs$SampleGroup, i, j) )
# between-group sum of squares:
ssa <- sst - ssw
msg <- c(paste("Adonis SumOfSqs for SampleGroup:", result$aov.tab$SumsOfSqs[1]),
         paste("Manual SumOfSqs for SampleGroup:", ssa))
cat(msg, sep="\n")

# Now, F-statistic:
# This is the F value for the samples within the groups we actually have.
# Permutations will follow.
f <- (ssa / (num_groups - 1)) / (ssw / (N*num_groups - num_groups))
msg <- c(paste("Adonis F.model for SampleGroup:", result$aov.tab$F.Model[1]),
         paste("Manual F.model for SampleGroup:", f))
cat(msg, sep="\n")

# And P-value:
# Here we need the permutations to shuffle items between groups.
# This does not currently work as it should

make_fp <- function(d, perms, g, grps, sst) {
  apply(perms, 1, function(perm) {
    a <- grps[perm]
    ssw2 <- g*diaggr(d, function(i, j) d[i, j]^2 * eps(a, i, j) )
    ssa2 <- sst - ssw2
    (ssa2 / (g - 1)) / (ssw2 / (N*g - g))
  })
}

f2 <- make_fp(d, perms, num_groups, attrs$SampleGroup, sst)

p <- (sum(f2 >= f) + 1) / (length(f2) + 1)
msg <- c(paste("Adonis p-value for SampleGroup:", result$aov.tab$`Pr(>F)`[1]),
         paste("Manual p-value for SampleGroup:", p))
cat(msg, sep="\n")


# Troubleshooting ---------------------------------------------------------

N <- 4
num_groups <- 2

dists <- mock_dist_mat(N, sd=0.001)
attrs <- data.frame(SampleGroup = c("A", "A", "A", "A", "B", "B", "B", "B"))
rownames(attrs) <- attr(dists, "Labels")

perms <- matrix(c(1:(N*num_groups),
                c(8,7,6,5,4,2,3,1),
                c(8,6,3,4,2,1,5,7)),
                nrow=3, byrow = T)

result <- adonis(dists ~ SampleGroup, data = attrs, permutations = perms)


d <- as.matrix(dists)
sst <- diaggr(d, function(i, j) d[i, j]^2 )
ssw <- num_groups*diaggr(d, function(i, j) d[i, j]^2 * eps(attrs$SampleGroup, i, j) )
ssa <- sst - ssw
f <- (ssa / (num_groups - 1)) / (ssw / (N*num_groups - num_groups))


f2 <- make_fp(d, perms, num_groups, attrs$SampleGroup, sst)

f2
as.vector(result$f.perms)