# Differential expression methods for bulk RNA-seq, benchmarked on simulated
# counts with known truth.
#
# Simulation: 6,000 genes with realistic mean expression and a mean-dependent
# negative binomial dispersion; 10% of genes differ between two groups with
# absolute log2 fold changes between 0.5 and 2; sample-specific library sizes.
# Scenarios vary the number of replicates per group (3, 5, 10), add a
# complete null (no true differences), and add outlier counts in one sample.
#
# Methods: DESeq2 (Wald test, default independent filtering and Cook's
# distance outlier handling), edgeR quasi-likelihood F-test, limma-voom,
# Welch t-test on log CPM, and Wilcoxon rank-sum on CPM.
#
# Usage: Rscript R/benchmark.R [replicates]

suppressPackageStartupMessages({
  library(DESeq2); library(edgeR); library(limma); library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(args)) as.integer(args[1]) else 8
G <- 6000; P_DE <- 0.10
out_tab <- "results/tables"; out_fig <- "results/figures"

simulate <- function(n, seed, de = TRUE, outlier = FALSE) {
  set.seed(seed)
  mu <- exp(rnorm(G, log(150), 1.7))                          # mean counts across genes
  disp <- (0.03 + 1.5 / mu) * exp(rnorm(G, 0, 0.5))           # dispersion falls with expression
  is_de <- rep(FALSE, G)
  lfc <- rep(0, G)
  if (de) {
    is_de[sample(G, P_DE * G)] <- TRUE
    lfc[is_de] <- sample(c(-1, 1), sum(is_de), TRUE) * runif(sum(is_de), 0.5, 2)
  }
  group <- factor(rep(c("A", "B"), each = n))
  sf <- exp(rnorm(2 * n, 0, 0.25))
  m <- outer(mu, sf) * ifelse(outer(rep(1, G), group == "B"), 2^lfc, 1)
  counts <- matrix(rnbinom(length(m), mu = m, size = 1 / disp), nrow = G)
  if (outlier) {                                              # one sample, 1% of genes inflated tenfold
    k <- sample(G, G / 100)
    counts[k, 1] <- counts[k, 1] * 10 + 20
  }
  rownames(counts) <- paste0("g", seq_len(G))
  list(counts = counts, group = group, truth = is_de)
}

run_methods <- function(s) {
  cts <- s$counts; grp <- s$group; p <- list()
  dds <- DESeqDataSetFromMatrix(cts, data.frame(group = grp), ~group)
  dds <- DESeq(dds, quiet = TRUE)
  p$`DESeq2` <- results(dds, contrast = c("group", "B", "A"))$padj

  y <- DGEList(cts, group = grp)
  keep <- filterByExpr(y)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  design <- model.matrix(~grp)
  fit <- glmQLFit(estimateDisp(y, design), design)
  q <- p.adjust(glmQLFTest(fit, coef = 2)$table$PValue, "BH")
  p$`edgeR QL` <- setNames(rep(NA, G), rownames(cts)); p$`edgeR QL`[keep] <- q

  v <- voom(y, design)
  q <- p.adjust(eBayes(lmFit(v, design))$p.value[, 2], "BH")
  p$`limma-voom` <- setNames(rep(NA, G), rownames(cts)); p$`limma-voom`[keep] <- q

  lcpm <- cpm(y, log = TRUE)
  tt <- apply(lcpm, 1, function(x) t.test(x[grp == "B"], x[grp == "A"])$p.value)
  p$`t-test on log CPM` <- setNames(rep(NA, G), rownames(cts)); p$`t-test on log CPM`[keep] <- p.adjust(tt, "BH")

  c2 <- cpm(y)
  wt <- apply(c2, 1, function(x) suppressWarnings(wilcox.test(x[grp == "B"], x[grp == "A"])$p.value))
  p$`Wilcoxon on CPM` <- setNames(rep(NA, G), rownames(cts)); p$`Wilcoxon on CPM`[keep] <- p.adjust(wt, "BH")
  p
}

score <- function(p, truth) {
  do.call(rbind, lapply(names(p), function(m) {
    called <- !is.na(p[[m]]) & p[[m]] < 0.05
    data.frame(method = m, called = sum(called), false = sum(called & !truth),
               fdp = if (sum(called)) sum(called & !truth) / sum(called) else 0,
               power = if (any(truth)) sum(called & truth) / sum(truth) else NA)
  }))
}

scenarios <- rbind(
  data.frame(scenario = "10% DE", n = c(3, 5, 10), de = TRUE, outlier = FALSE),
  data.frame(scenario = "no DE (null)", n = c(3, 5, 10), de = FALSE, outlier = FALSE),
  data.frame(scenario = "10% DE, outlier sample", n = 5, de = TRUE, outlier = TRUE))

res <- list()
for (i in seq_len(nrow(scenarios))) {
  sc <- scenarios[i, ]
  for (r in seq_len(REPS)) {
    s <- simulate(sc$n, seed = 1000 * i + r, de = sc$de, outlier = sc$outlier)
    out <- score(run_methods(s), s$truth)
    res[[length(res) + 1]] <- cbind(scenario = sc$scenario, n_per_group = sc$n, rep = r, out)
  }
  message(sc$scenario, ", n = ", sc$n, " done")
}
res <- do.call(rbind, res)
write.csv(res, file.path(out_tab, "de_benchmark_all_runs.csv"), row.names = FALSE)
summ <- aggregate(cbind(called, false, fdp, power) ~ scenario + n_per_group + method, data = res,
                  FUN = mean, na.action = na.pass)
summ[, c("called", "false")] <- round(summ[, c("called", "false")], 1)
summ[, c("fdp", "power")] <- round(summ[, c("fdp", "power")], 3)
summ <- summ[order(summ$scenario, summ$n_per_group, summ$method), ]
write.csv(summ, file.path(out_tab, "de_benchmark_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)

# Under the complete null every call is false; count the runs with any.
nul <- subset(res, scenario == "no DE (null)")
null_tab <- do.call(rbind, lapply(split(nul, list(nul$n_per_group, nul$method), drop = TRUE), function(d)
  data.frame(n_per_group = d$n_per_group[1], method = d$method[1], runs = nrow(d),
             runs_with_any_false_call = sum(d$false > 0), mean_false_calls = mean(d$false),
             max_false_calls = max(d$false))))
null_tab <- null_tab[order(null_tab$n_per_group, null_tab$method), ]
write.csv(null_tab, file.path(out_tab, "null_scenario.csv"), row.names = FALSE)

d <- subset(summ, scenario == "10% DE")
long <- rbind(data.frame(d[, c("n_per_group", "method")], metric = "False discovery proportion (target 0.05)", value = d$fdp),
              data.frame(d[, c("n_per_group", "method")], metric = "Power (share of true DE genes found)", value = d$power))
g <- ggplot(long, aes(n_per_group, value, colour = method)) + geom_line() + geom_point() +
  facet_wrap(~metric, scales = "free_y") + scale_x_continuous(breaks = c(3, 5, 10)) +
  geom_hline(data = data.frame(metric = "False discovery proportion (target 0.05)", y = 0.05),
             aes(yintercept = y), linetype = "dotted") +
  labs(x = "Replicates per group", y = NULL, colour = NULL,
       title = "Differential expression methods on simulated counts with 10% true DE genes") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out_fig, "fdp_power_by_sample_size.png"), g, width = 10, height = 4.5, dpi = 150)

o <- subset(summ, scenario %in% c("10% DE", "10% DE, outlier sample") & n_per_group == 5)
g2 <- ggplot(o, aes(method, fdp, fill = scenario)) + geom_col(position = "dodge") +
  geom_hline(yintercept = 0.05, linetype = "dotted") +
  labs(x = NULL, y = "False discovery proportion", fill = NULL,
       title = "Effect of outlier counts in one of ten samples (5 per group)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(out_fig, "outlier_effect.png"), g2, width = 8, height = 4.5, dpi = 150)
