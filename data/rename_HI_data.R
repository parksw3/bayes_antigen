library(seqinr)

rr1 <- read.fasta("../data/gisaid_epiflu_sequence.fasta", seqtype="AA", as.string=TRUE)
rr2 <- read.fasta("../data/ProteinFastaResults.fasta", seqtype="AA", as.string=TRUE)

seq1 <- sapply(rr1, "[[", 1)
names(seq1) <- gsub(" \\|.*", "", gsub(">", "", unname(sapply(rr1, function(x) attr(x, "Annot")))))
seq1_accession <- gsub(".*\\| ", "", gsub(">", "", unname(sapply(rr1, function(x) attr(x, "Annot")))))

seq2 <- sapply(rr2, "[[", 1)
name2 <- unname(gsub("Protein.*", "", gsub(".*Strain Name:", "", sapply(rr2, function(x) attr(x, "Annot")))))
names(seq2) <- substr(name2, 1, nchar(name2)-1) 
seq2_accession <-  gsub("\\|.*", "", gsub(">gb:", "", sapply(rr2, function(x) attr(x, "Annot"))))

H1N1_data <- read.table("../data/H1N1_HI_data.tsv", sep = '\t', header = TRUE)

seq_data <- read.table("../data/H1N1_seq_data.tsv", sep = '\t', header = TRUE)

name_data <- data.frame(
	accession=unname(c(seq1_accession, seq2_accession)),
	new_names=c(names(seq1), names(seq2))
)

new_seq_data <- merge(seq_data, name_data)

## which missing
seq_data[!(seq_data$accession %in% name_data$accession),]

H1N1_data$virusStrain <- new_seq_data$new_names[match(H1N1_data$virusStrain, new_seq_data$strain)]
H1N1_data$serumStrain <- new_seq_data$new_names[match(H1N1_data$serumStrain, new_seq_data$strain)]

H1N1_new_data <- H1N1_data[!is.na(H1N1_data$virusStrain) & !is.na(H1N1_data$serumStrain),]

save("H1N1_new_data", "new_seq_data", file="rename_HI_data.rda")
