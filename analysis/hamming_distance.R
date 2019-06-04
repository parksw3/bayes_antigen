library(seqinr)
library(DECIPHER)

rr1 <- read.fasta("../data/gisaid_epiflu_sequence.fasta", seqtype="AA", as.string=TRUE)
rr2 <- read.fasta("../data/ProteinFastaResults.fasta", seqtype="AA", as.string=TRUE)

seq1 <- sapply(rr1, "[[", 1)
names(seq1) <- gsub(" \\|.*", "", gsub(">", "", unname(sapply(rr1, function(x) attr(x, "Annot")))))

seq2 <- sapply(rr2, "[[", 1)
name2 <- unname(gsub("Protein.*", "", gsub(".*Strain Name:", "", sapply(rr2, function(x) attr(x, "Annot")))))
names(seq2) <- substr(name2, 1, nchar(name2)-1) 

all_seq <- AlignSeqs(AAStringSet(c(seq1, seq2)))

nseq <- length(all_seq)

sa_mat <- sb_mat <- ca1_mat <- ca2_mat <- cb_mat <- all_mat <- matrix(NA, nseq, nseq)

sa_region <- c(141:142, 170:174, 176:181)
sb_region <- c(201:212)
ca1_region <- c(183:187, 220:222, 252:254)
ca2_region <- c(154:159, 238:239)
cb_region <- c(87:92)

all_region <- c(sa_region, sb_region, ca1_region, ca2_region, cb_region)
## checking whether anything is missing... 
which(unname(sapply(all_seq, function(x) any(strsplit(as.character(x[all_region]), "")[[1]]=="-"))))

## missing cb_region...
strsplit(as.character(all_seq[60]), "")[[1]][all_region]

distfun <- function(sp1, sp2, region) {
	mm <- sp1[region]!=sp2[region]
	
	complete <- ((sp1[region]!="-") & (sp2[region]!="-"))
	
	sum(mm[complete])/length(mm[complete])
}

for (i in 1:nseq) {
	print(i)
	for (j in i:nseq) {
		ss1 <- as.character(all_seq[[i]])
		ss2 <- as.character(all_seq[[j]])
		
		sp1 <- strsplit(ss1, "")[[1]]
		sp2 <- strsplit(ss2, "")[[1]]
		
		sa_mat[i,j] <- sa_mat[j,i] <- distfun(sp1, sp2, sa_region)
		sb_mat[i,j] <- sb_mat[j,i] <- distfun(sp1, sp2, sb_region)
		ca1_mat[i,j] <- ca1_mat[j,i] <- distfun(sp1, sp2, ca1_region)
		ca2_mat[i,j] <- ca2_mat[j,i] <- distfun(sp1, sp2, ca2_region)
		cb_mat[i,j] <- cb_mat[j,i] <- distfun(sp1, sp2, cb_region)
	}
}

rownames(sa_mat) <- colnames(sa_mat) <-
	rownames(sb_mat) <- colnames(sb_mat) <-
	rownames(ca1_mat) <- colnames(ca1_mat) <-
	rownames(ca2_mat) <- colnames(ca2_mat) <-
	rownames(cb_mat) <- colnames(cb_mat) <-
	names(all_seq)

save("sa_mat", "sb_mat", "ca1_mat", "ca2_mat", "cb_mat", "all_seq", file="hamming_distance.rda")
