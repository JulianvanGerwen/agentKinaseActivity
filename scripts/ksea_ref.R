# KSEAapp reference for kinact.methods.ksea. Usage: Rscript ksea_ref.R lib.csv px_dir out_dir
# lib.csv: kinase, acc, rsd ; px_dir: one CSV per condition with columns Protein,Gene,Peptide,Residue.Both,p,FC
suppressMessages(library(KSEAapp))
a <- commandArgs(TRUE)
data(KSData)
lib <- read.csv(a[1], stringsAsFactors = FALSE)
ks <- KSData[rep(1, nrow(lib)), ]
ks$KINASE <- lib$kinase; ks$GENE <- lib$kinase; ks$KIN_ACC_ID <- "X"; ks$SUB_ACC_ID <- lib$acc
ks$SUB_GENE <- lib$acc; ks$SUB_MOD_RSD <- lib$rsd; ks$Source <- "PhosphoSitePlus"; ks$networkin_score <- Inf
rownames(ks) <- NULL
for (f in list.files(a[2], pattern = "csv$", full.names = TRUE)) {
  px <- read.csv(f, stringsAsFactors = FALSE)
  sc <- KSEA.Scores(ks, px, NetworKIN = FALSE, NetworKIN.cutoff = 5)
  write.csv(sc, file.path(a[3], basename(f)), row.names = FALSE)
}
