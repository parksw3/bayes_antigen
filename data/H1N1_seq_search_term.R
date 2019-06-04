library(dplyr)

rr <- read.table("H1N1_seq_data.tsv",
				 header=TRUE)

rr %>% 
	filter(database=="GISAID") %>%
	select(accession) %>%
	unlist %>%
	as.character %>%
	gsub(pattern="EPI_ISL_", replacement="") %>% 
	paste(collapse=" ")

rr %>% 
	filter(database!="GISAID") %>%
	select(strain) %>%
	unlist %>%
	as.character %>%
	paste(collapse=", ")

rr %>% 
	filter(database!="GISAID") %>%
	select(strain, accession) 

