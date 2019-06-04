library(dplyr)
library(bayesplot)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

wquant <- function (x, weights, probs = c(0.025, 0.975)) {
	idx <- order(x)
	x <- x[idx]
	weights <- weights[idx]
	w <- cumsum(weights)/sum(weights)
	rval <- approx(w,x,probs,rule=1)
	rval$y
}

load("../analysis/h1n1_analysis_relax2.rda")

ee_h1n1 <- rstan::extract(H1N1_fit)

h1n1_virus1 <- data_frame(
	x=apply(ee_h1n1$X[1:500,1,], 2, median),
	y=apply(ee_h1n1$X[1:500,2,], 2, median),
	virus=levels(H1N1_mod$virusStrain),
	virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)]
)

h1n1_serum1 <- data_frame(
	x=apply(ee_h1n1$Y[1:500,1,], 2, median),
	y=apply(ee_h1n1$Y[1:500,2,], 2, median),
	serum=levels(H1N1_mod$serumStrain),
	serumYear=H1N1_mod$serumYear[match(serum, H1N1_mod$serumStrain)]
)

h1n1_virus2 <- data_frame(
	x=apply(ee_h1n1$X[1:500+500,1,], 2, median),
	y=apply(ee_h1n1$X[1:500+500,2,], 2, median),
	virus=levels(H1N1_mod$virusStrain),
	virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)]
)

h1n1_serum2 <- data_frame(
	x=apply(ee_h1n1$Y[1:500+500,1,], 2, median),
	y=apply(ee_h1n1$Y[1:500+500,2,], 2, median),
	serum=levels(H1N1_mod$serumStrain),
	serumYear=H1N1_mod$serumYear[match(serum, H1N1_mod$serumStrain)]
)

h1n1_virus3 <- data_frame(
	x=apply(ee_h1n1$X[1:500+1000,1,], 2, median),
	y=apply(ee_h1n1$X[1:500+1000,2,], 2, median),
	virus=levels(H1N1_mod$virusStrain),
	virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)]
)

h1n1_serum3 <- data_frame(
	x=apply(ee_h1n1$Y[1:500+1000,1,], 2, median),
	y=apply(ee_h1n1$Y[1:500+1000,2,], 2, median),
	serum=levels(H1N1_mod$serumStrain),
	serumYear=H1N1_mod$serumYear[match(serum, H1N1_mod$serumStrain)]
)

h1n1_virus4 <- data_frame(
	x=apply(ee_h1n1$X[1:500+1500,1,], 2, median),
	y=apply(ee_h1n1$X[1:500+1500,2,], 2, median),
	virus=levels(H1N1_mod$virusStrain),
	virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)]
) 

h1n1_serum4 <- data_frame(
	x=apply(ee_h1n1$Y[1:500+1500,1,], 2, median),
	y=apply(ee_h1n1$Y[1:500+1500,2,], 2, median),
	serum=levels(H1N1_mod$serumStrain),
	serumYear=H1N1_mod$serumYear[match(serum, H1N1_mod$serumStrain)]
)

g1 <- ggplot(h1n1_virus1) +
	geom_point(aes(x, y, col=virusYear)) +
	xlab("Antigenic dimension 1") +
	ylab("Antigenic dimension 2") +
	scale_color_gradient(low="red", high="black") +
	ggtitle("Chain 1") +
	theme(
		legend.position = "none"
	)

g2 <- (g1 %+% h1n1_virus2) +
	ggtitle("Chain 2")

g3 <- (g1 %+% h1n1_virus3) +
	ggtitle("Chain 3")

g4 <- (g1 %+% h1n1_virus4) +
	ggtitle("Chain 4")

gmap <- arrangeGrob(
	g1 + geom_point(data=h1n1_serum1, aes(x, y), shape=2, col="blue"), 
	g2 + geom_point(data=h1n1_serum2, aes(x, y), shape=2, col="blue"), 
	g3 + geom_point(data=h1n1_serum3, aes(x, y), shape=2, col="blue"),
	g4 + geom_point(data=h1n1_serum4, aes(x, y), shape=2, col="blue"), 
	nrow=2)

ggsave("h1n1_map_relax2.pdf", gmap, width=8, height=6)

H1N1_rhat <- bayesplot::rhat(H1N1_fit)
hist(H1N1_rhat[xg <- grepl("X", names(H1N1_rhat))])
hist(H1N1_rhat[yg <- grepl("Y", names(H1N1_rhat))])
hist(H1N1_rhat[!(xg | yg)])

distance_data <- lapply(1:2000, function(n) {
	X <- ee_h1n1$X[n,,]
	ref <- X[,112]
	
	data_frame(
		virus=levels(H1N1_mod$virusStrain),
		virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)],
		dist=sqrt(colSums((X - ref)^2))
	)
}) %>%
	bind_rows(.id="iter")

distance_summary <- distance_data %>%
	mutate(chain=ceiling(as.numeric(iter)/500)) %>%
	mutate(weight=ifelse(chain %in% c(1, 3), 1/2, 1)) %>%
	group_by(virusYear, virus) %>%
	summarize(
		median=wquant(dist, weight, 0.5),
		lwr=wquant(dist, weight, 0.025),
		upr=wquant(dist, weight, 0.975)
	)
	
gdist <- ggplot(distance_summary) +
	geom_errorbar(aes(virusYear, ymin=lwr, ymax=upr, group=virus), 
				  width=0, position = position_dodge(width=0.4), col="grey") +
	geom_point(aes(virusYear, median, group=virus), position = position_dodge(width=0.4)) +
	geom_smooth(aes(virusYear, median), col="red", se=FALSE) +
	scale_x_continuous("Year") +
	scale_y_continuous("Antigenic distance") +
	theme(
		strip.background = element_blank(),
		panel.spacing = grid::unit(0, "cm")
	)

ggsave("h1n1_distance_relax2.pdf", gdist, width=8, height=6)

w <- c(rep(1, 500), rep(1, 500), rep(1, 500), rep(1, 500))

h1n1_virus_effect <- data_frame(
	median=apply(ee_h1n1$re_J, 2, wquant, weights=w, 0.5),
	lwr=apply(ee_h1n1$re_J, 2, wquant, weights=w, 0.025),
	upr=apply(ee_h1n1$re_J, 2, wquant, weights=w, 0.975),
	virus=levels(H1N1_mod$virusStrain),
	year=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)],
	key="A/H1N1"
)

h1n1_serum_effect <- data_frame(
	median=apply(ee_h1n1$re_A, 2, wquant, weights=w, 0.5),
	lwr=apply(ee_h1n1$re_A, 2, wquant, weights=w, 0.025),
	upr=apply(ee_h1n1$re_A, 2, wquant, weights=w, 0.975),
	serum=levels(H1N1_mod$serumStrain),
	year=H1N1_mod$serumYear[match(serum, H1N1_mod$serumStrain)],
	key="A/H1N1"
)

gvirus <- ggplot(h1n1_virus_effect) +
	geom_hline(yintercept=0, lty=2) +
	geom_boxplot(aes(year, median, group=interaction(year)), alpha=0.1, varwidth=TRUE) +
	geom_smooth(aes(year, median), span=0.2, se=FALSE, lwd=1, n=300, col="black") +
	scale_y_continuous("Effects of virus strain") +
	scale_x_continuous("Year") +
	facet_grid(key~., scale="free_y") +
	theme(
		panel.grid = element_blank(),
		panel.spacing = grid::unit(0, "cm"),
		legend.position = "none",
		strip.background = element_blank()
	)

gserum <- gvirus %+% h1n1_serum_effect +
	scale_y_continuous("Effects of serum strain")

gtot <- arrangeGrob(gvirus, gserum, ncol=1)

ggsave("h1n1_effect_relax2.pdf", gtot, width=6, height=4)
