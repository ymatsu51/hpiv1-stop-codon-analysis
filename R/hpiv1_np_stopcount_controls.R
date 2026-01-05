## ---- packages ----
pkgs_cran <- c("rentrez", "Biostrings")
for (p in pkgs_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(rentrez)
library(Biostrings)

## ---- parameters ----
term   <- '"Human respirovirus 1"[Organism] AND "complete genome"[Title]'
retmax <- 5000

## Insert between "from-the-end 758 | 757" (i.e., after the base that is 758th from the end)
## For each sequence length L (1-based), this corresponds to pos_each = L - 758 + 1
K_FROM_END_AFTER <- 758L
INSERT_NT <- "G"

out_all_cds_fa <- "HPIV1_all_CDS.fasta"
out_stats      <- "HPIV1_NP_edited_stop_stats.tsv"

## ---- helpers (S4-free) ----
dna_len <- function(dna) nchar(as.character(dna))
aa_len  <- function(aa)  nchar(as.character(aa))

sanitize_dna <- function(dna) {
  s <- toupper(as.character(dna))
  s <- gsub("U", "T", s)         # precautionary conversion
  s <- gsub("[^ACGTN]", "N", s)  # convert non-ACGTN characters to N (gaps, ?, spaces, etc.)
  Biostrings::DNAString(s)
}

insert_nt_after <- function(dna, pos, nt = "G") {
  s <- as.character(dna)
  L <- nchar(s)
  if (pos < 0L || pos >= L) stop("pos must satisfy 0 <= pos < length(dna). Got pos=", pos, " L=", L)
  Biostrings::DNAString(paste0(substr(s, 1L, pos), nt, substr(s, pos + 1L, L)))
}

stop_positions <- function(aa) {
  m <- gregexpr("\\*", as.character(aa))[[1]]
  if (length(m) == 1L && m[1] == -1L) integer(0) else as.integer(m)
}

## Identify NP (N gene) CDS entries from FASTA headers (robust to header variation)
## - Prioritize gene=N annotation (N is short, so gene context is important)
## - Also capture product names such as nucleoprotein / nucleocapsid (protein)
is_np_cds_header <- function(h) {
  h2 <- tolower(h)

  geneN <- grepl("gene[= :]+n([^a-z0-9]|$)", h2) ||
           grepl("\\bgene\\b[^\\n]*\\bn\\b", h2)

  prodN <- grepl("nucleoprotein", h2) ||
           grepl("nucleocapsid", h2) ||
           grepl("\\bn protein\\b", h2) ||
           grepl("nucleocapsid protein", h2)

  geneN || prodN
}

## Calculate the insertion position (1-based) for each sequence:
## pos_each = L - K_FROM_END_AFTER + 1
## insert_nt_after() inserts between pos|pos+1, so pos_each can be used directly
calc_pos_each <- function(L, k_from_end_after) {
  as.integer(L - k_from_end_after + 1L)
}

## =============================================================================
## Step A) Retrieve CDS FASTA from NCBI and extract NP (N) CDS
## =============================================================================
srch <- rentrez::entrez_search(db = "nuccore", term = term, retmax = retmax)
cat("Hits:", srch$count, "Fetched IDs:", length(srch$ids), "\n")
if (length(srch$ids) == 0) stop("No records found for query: ", term)

cds_fasta_txt <- rentrez::entrez_fetch(
  db = "nuccore",
  id = srch$ids,
  rettype = "fasta_cds_na",
  retmode = "text"
)
writeLines(cds_fasta_txt, out_all_cds_fa)

cds_all <- Biostrings::readDNAStringSet(out_all_cds_fa)
cat("Total CDS records fetched:", length(cds_all), "\n")
if (length(cds_all) == 0) stop("No CDS FASTA entries were fetched.")

keepN <- vapply(names(cds_all), is_np_cds_header, logical(1))
np_nt <- cds_all[keepN]
cat("NP(N)-like CDS records:", length(np_nt), "\n")
if (length(np_nt) == 0) stop("No NP-like CDS found. Inspect FASTA headers and adjust is_np_cds_header().")

## Retain only sequences long enough to define the insertion site
## (at least K_FROM_END_AFTER nucleotides required)
L_nt <- vapply(as.list(np_nt), dna_len, integer(1))
np_nt_ok <- np_nt[L_nt >= K_FROM_END_AFTER]
cat("Usable NP CDS (len >= ", K_FROM_END_AFTER, " nt): ", length(np_nt_ok), "\n", sep = "")
if (length(np_nt_ok) == 0) stop("No NP CDS sequences are long enough for insertion at from-end position.")

## =============================================================================
## Step B) Sanitize + single-nucleotide insertion (from-end 758|757)
##          + sanitize + fuzzy translation
## =============================================================================

## 1) Sanitize original CDS
np_nt_ok_s <- Biostrings::DNAStringSet(lapply(as.list(np_nt_ok), sanitize_dna))
names(np_nt_ok_s) <- names(np_nt_ok)

## 2) Calculate pos_each for each sequence
L_ok <- vapply(as.list(np_nt_ok_s), dna_len, integer(1))
pos_each <- vapply(L_ok, calc_pos_each, integer(1), k_from_end_after = K_FROM_END_AFTER)

## Validate insertion positions: 0 <= pos < L
bad <- which(pos_each < 0L | pos_each >= L_ok)
if (length(bad) > 0) {
  cat("Found invalid insertion positions for", length(bad), "sequences.\n")
  cat("Example headers:\n")
  print(head(names(np_nt_ok_s)[bad], 5))
  stop("Invalid insertion position(s) detected. Please inspect sequences/lengths.")
}

## 3) Insert one nucleotide (sequence-specific positions via mapply)
np_nt_edit_list <- mapply(
  FUN = function(dna, pos) insert_nt_after(dna, pos = pos, nt = INSERT_NT),
  dna = as.list(np_nt_ok_s),
  pos = pos_each,
  SIMPLIFY = FALSE
)
np_nt_edit <- Biostrings::DNAStringSet(np_nt_edit_list)

## Append insertion position information to FASTA headers
## - From the 5' end: pos|pos+1
## - From the 3' end target (e.g., fromEnd=758|757)
names(np_nt_edit) <- paste0(
  names(np_nt_ok_s),
  "|Gins_pos=", pos_each, "|", pos_each + 1L,
  "|fromEnd=", K_FROM_END_AFTER, "|", K_FROM_END_AFTER - 1L
)

## 4) Sanitize edited sequences as a precaution
np_nt_edit_s <- Biostrings::DNAStringSet(lapply(as.list(np_nt_edit), sanitize_dna))
names(np_nt_edit_s) <- names(np_nt_edit)

## 5) Fuzzy translation: ambiguous codons translated as X (do not terminate)
aa_native <- Biostrings::translate(np_nt_ok_s,  if.fuzzy.codon = "X")
aa_edited <- Biostrings::translate(np_nt_edit_s, if.fuzzy.codon = "X")
names(aa_native) <- names(np_nt_ok_s)
names(aa_edited) <- names(np_nt_edit_s)

## =============================================================================
## Step C) Generate stop-codon structure TSV (main output)
## =============================================================================
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
  cds_len_nt_native = as.integer(L_ok),
  cds_len_nt_edited = as.integer(vapply(as.list(np_nt_edit_s), dna_len, integer(1))),
  insert_after_pos_1based = as.integer(pos_each),
  insert_between_pos_1based = paste0(pos_each, "|", pos_each + 1L),
  from_end_target = paste0(K_FROM_END_AFTER, "|", K_FROM_END_AFTER - 1L),
  aa_len_native = as.integer(aaL_native),
  aa_len_edited = as.integer(aaL_edited),
  stop_count = as.integer(stop_n),
  stop_density = as.numeric(stop_density),
  max_stop_free_run = as.integer(max_stop_free),
  stop_positions_aa_1based = stop_positions_str,
  stringsAsFactors = FALSE
)

utils::write.table(stats, out_stats, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Wrote TSV:", out_stats, "\n",
    "N sequences:", nrow(stats), "\n", sep = "")
