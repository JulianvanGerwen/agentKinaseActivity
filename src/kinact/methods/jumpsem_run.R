# Run JUMPsem 1.1.0 natively on one study. Usage: Rscript jumpsem_run.R input.csv db.csv out_activity.tsv
# input.csv: Substrate_ACC, Substrate_GN, sites, Peptide, S1..Sn (log2FC, complete cases); db.csv: PSP-format table.
suppressMessages({library(JUMPsem)})
args <- commandArgs(TRUE)
inp <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
db <- read.csv(args[2], check.names = FALSE, stringsAsFactors = FALSE)
od <- dirname(args[3])
t0 <- Sys.time()
res <- suppressWarnings(suppressMessages(
  JUMPsem(input = inp, datatype = "psp", organism = "human", database = db, input.log2.norm = FALSE,
          output.folder = od)))
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
if (inherits(res, "try-error") || is.null(res$Activity)) {
  writeLines(paste0("#error\t", el), args[3])
} else {
  write.table(res$Activity, args[3], sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(as.character(el), paste0(args[3], ".time"))
}
