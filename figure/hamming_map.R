library(dplyr)
library(ggplot2); theme_set(theme_bw())

load("../analysis/hamming_distance.rda")

sa_map <- cmdscale(sa_mat * 20)
sb_map <- cmdscale(sb_mat * 20)
ca1_map <- cmdscale(ca1_mat * 20)
ca2_map <- cmdscale(ca2_mat * 20)
cb_map <- cmdscale(cb_mat[-60, -60] * 20)

allmap <- list(Sa=sa_map, Sb=sb_map, Ca1=ca1_map, Ca2=ca2_map, Cb=cb_map) %>%
	lapply(function(x){
		nn <- rownames(x)
		year <- gsub(".*/", "", nn)
		year <- gsub(" .*", "", year)
		year <- as.numeric(unname(sapply(year, function(x){ 
			gsub(" ", "", ifelse(length(strsplit(x, "")[[1]])==2, paste("19", x, collapse=""), x))
		})))
		
		year[year==1905] <- 2005
		
		dd <- as.data.frame(x)
		dd$name <- nn
		dd$year <- year
		dd
	}) %>%
	bind_rows(.id="region")

allmap$id <- 1:nrow(allmap)

g1 <- ggplot(allmap) +
	geom_point(aes(V1, V2, col=year, group=id), position=position_jitter(width=0.1, height=0.1)) +
	scale_color_gradient(low="red", high="blue") +
	xlab("Dimension 1") +
	ylab("Dimension 2") +
	facet_wrap(~region)

cb_mat[is.na(cb_mat)] <- 0

combined_map <- cmdscale((sa_mat + sb_mat + ca1_mat + ca2_mat + cb_mat)/5 * 20)
nn <- rownames(combined_map)
combined_map <- as.data.frame(combined_map)

combined_map$year <- allmap$year[match(nn, allmap$name)]
combined_map$name <- nn

g2 <- ggplot(combined_map) +
	geom_point(aes(V1, V2, col=year)) +
	scale_color_gradient(low="red", high="blue") +
	xlab("Dimension 1") +
	ylab("Dimension 2")

ggsave("individual_map.pdf", g1, width=8, height=6)
ggsave("average_map.pdf", g2, width=8, height=6)
save("combined_map", file="hamming_map.rda")
