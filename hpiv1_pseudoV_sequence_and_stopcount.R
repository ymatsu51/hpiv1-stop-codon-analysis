## -> search complete genomes
## -> fetch CDS FASTA (fasta_cds_na)
## -> keep P CDS
## -> insert 1G between 949|950 (after position 949)
## -> translate full length (keep stops)
## -> enumerate ALL stop (*) positions and summarize
##
## Outputs:
##   - HPIV1_all_CDS.fasta
##   - HPIV1_P_CDS_native.fasta
##   - HPIV1_P_CDS_edited_Gins949_950.fasta
##   - HPIV1_P_native_AA.fasta
##   - HPIV1_P_edited_AA_full_with_stops.fasta
##   - HPIV1_P_edited_stop_stats.tsv

## ---- packages ----
pkgs_cran <- c("rentrez", "Biostrings")
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(rentrez)
library(Biostrings)

## ---- parameters ----
term <- '"Human respirovirus 1"[Organism] AND "complete genome"[Title]'
retmax <- 5000

INS_POS <- 949L  # insert AFTER this position: between 949|950 (1-based; ATG A=1)

out_all_cds_fa <- "HPIV1_all_CDS.fasta"

out_nt_native <- "HPIV1_P_CDS_native.fasta"
out_nt_edited <- "HPIV1_P_CDS_edited_Gins949_950.fasta"
out_aa_native <- "HPIV1_P_native_AA.fasta"
out_aa_edited <- "HPIV1_P_edited_AA_full_with_stops.fasta"
out_stats     <- "HPIV1_P_edited_stop_stats.tsv"

## ---- helpers (S4-free length/stop handling) ----
dna_len <- function(dna) nchar(as.character(dna))
aa_len  <- function(aa)  nchar(as.character(aa))

insert_nt_after <- function(dna, pos, nt = "G") {
  # pos: insert AFTER pos, i.e., between pos|pos+1
  s <- as.character(dna)
  L <- nchar(s)
  if (pos < 0L || pos >= L) stop("pos must satisfy 0 <= pos < length(dna). Got pos=", pos, " L=", L)
  Biostrings::DNAString(paste0(substr(s, 1L, pos), nt, substr(s, pos + 1L, L)))
}

stop_positions <- function(aa) {
  # 1-based positions of '*' in an AAString
  m <- gregexpr("\\*", as.character(aa))[[1]]
  if (length(m) == 1L && m[1] == -1L) integer(0) else as.integer(m)
}

## Tolerant P-CDS header filter (NCBI fasta_cds_na headers vary)
is_p_cds_header <- function(h) {
  h2 <- tolower(h)

  # gene annotation patterns (varies across NCBI headers)
  geneP <- grepl("gene[= :]+p([^a-z0-9]|$)", h2) ||
           grepl("\\bgene\\b[^\\n]*\\bp\\b", h2)

  # product/protein name patterns
  prodP <- grepl("phosphoprotein", h2) ||
           grepl("\\bp protein\\b", h2)

  geneP || prodP
}

## ---- 1) search complete genomes ----
srch <- rentrez::entrez_search(db = "nuccore", term = term, retmax = retmax)
cat("Hits:", srch$count, "Fetched IDs:", length(srch$ids), "\n")
if (length(srch$ids) == 0) stop("No records found for query: ", term)

## ---- 2) fetch CDS FASTA (nucleotide) from those genomes ----
cds_fasta_txt <- rentrez::entrez_fetch(
  db = "nuccore",
  id = srch$ids,
  rettype = "fasta_cds_na",
  retmode = "text"
)
writeLines(cds_fasta_txt, out_all_cds_fa)

## ---- 3) read & keep only P CDS ----
cds_all <- Biostrings::readDNAStringSet(out_all_cds_fa)
cat("Total CDS records fetched:", length(cds_all), "\n")
if (length(cds_all) == 0) stop("No CDS FASTA entries were fetched. Check query / NCBI availability / network.")

keepP <- vapply(names(cds_all), is_p_cds_header, logical(1))
p_nt <- cds_all[keepP]
cat("P-like CDS records:", length(p_nt), "\n")

if (length(p_nt) == 0) {
  cat("\nNo P-like CDS found with current header rules.\n",
      "Example headers (first 30):\n", sep = "")
  cat(paste(head(names(cds_all), 30), collapse = "\n"), "\n\n")
  stop("Adjust is_p_cds_header() to match your headers (e.g., gene=, product=).")
}

## Save native P CDS
Biostrings::writeXStringSet(p_nt, out_nt_native, width = 60)

## ---- 4) QC: sequences must be long enough for boundary 949|950 ----
L_nt <- vapply(as.list(p_nt), dna_len, integer(1))
keep_len <- L_nt >= (INS_POS + 1L)  # need at least 950 nt to have 949|950 boundary
p_nt_ok <- p_nt[keep_len]
cat("Usable P CDS (len >= 950 nt):", length(p_nt_ok), "\n")
if (length(p_nt_ok) == 0) stop("No P CDS sequences are long enough for insertion at 949|950.")

## ---- 5) make edited nt (+1G between 949|950) ----
p_nt_edit <- Biostrings::DNAStringSet(lapply(as.list(p_nt_ok), insert_nt_after, pos = INS_POS, nt = "G"))
names(p_nt_edit) <- paste0(names(p_nt_ok), "|Gins=", INS_POS, "|", INS_POS + 1L)

Biostrings::writeXStringSet(p_nt_edit, out_nt_edited, width = 60)

## ---- 6) translate full length (keep stops) ----
aa_native <- Biostrings::translate(p_nt_ok)    # may include terminal '*'
aa_edited <- Biostrings::translate(p_nt_edit)  # includes '*' wherever stops occur

names(aa_native) <- names(p_nt_ok)
names(aa_edited) <- names(p_nt_edit)

Biostrings::writeXStringSet(aa_native, out_aa_native, width = 60)
Biostrings::writeXStringSet(aa_edited, out_aa_edited, width = 60)

## ---- 7) enumerate ALL stop positions & summarize ----
aaL_native <- vapply(as.list(aa_native), aa_len, integer(1))
aaL_edited <- vapply(as.list(aa_edited), aa_len, integer(1))

stops_list <- lapply(as.list(aa_edited), stop_positions)

stop_n <- vapply(stops_list, length, integer(1))
stop_density <- stop_n / aaL_edited

max_stop_free <- mapply(function(pos, L) {
  if (length(pos) == 0) return(L)
  seg <- c(pos[1] - 1L, diff(pos) - 1L, L - pos[length(pos)])
  max(seg)
}, stops_list, aaL_edited)

stop_positions_str <- vapply(stops_list, function(pos) {
  if (length(pos) == 0) "" else paste(pos, collapse = ",")
}, character(1))

stats <- data.frame(
  header = names(aa_edited),
  cds_len_nt_native = as.integer(L_nt[keep_len]),
  cds_len_nt_edited = as.integer(L_nt[keep_len] + 1L),
  aa_len_native = as.integer(aaL_native),
  aa_len_edited = as.integer(aaL_edited),
  stop_count = as.integer(stop_n),
  stop_density = as.numeric(stop_density),
  max_stop_free_run = as.integer(max_stop_free),
  stop_positions_aa_1based = stop_positions_str,
  stringsAsFactors = FALSE
)

utils::write.table(stats, out_stats, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\n",
    "All CDS FASTA: ", out_all_cds_fa, "\n",
    "Native P nt:  ", out_nt_native, "\n",
    "Edited P nt:  ", out_nt_edited, "\n",
    "Native P aa:  ", out_aa_native, "\n",
    "Edited P aa:  ", out_aa_edited, "\n",
    "Stats:        ", out_stats, "\n", sep = "")
