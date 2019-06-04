library(dplyr)
library(bayesplot)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

load("../analysis/h1n1_analysis.rda")

ee_h1n1 <- rstan::extract(H1N1_fit)

h1n1_virus1 <- data_frame(
	x=apply(ee_h1n1$X[1:2000,1,], 2, median),
	y=apply(ee_h1n1$X[1:2000,2,], 2, median),
	virus=levels(H1N1_mod$virusStrain),
	virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)]
)

h1n1_serum1 <- data_frame(
	x=apply(ee_h1n1$Y[1:2000,1,], 2, median),
	y=apply(ee_h1n1$Y[1:2000,2,], 2, median),
	serum=levels(H1N1_mod$serumStrain),
	serumYear=H1N1_mod$serumYear[match(serum, H1N1_mod$serumStrain)]
)

h1n1_virus2 <- data_frame(
	x=apply(ee_h1n1$X[1:2000+2000,1,], 2, median),
	y=apply(ee_h1n1$X[1:2000+2000,2,], 2, median),
	virus=levels(H1N1_mod$virusStrain),
	virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)]
)

h1n1_serum2 <- data_frame(
	x=apply(ee_h1n1$Y[1:2000+2000,1,], 2, median),
	y=apply(ee_h1n1$Y[1:2000+2000,2,], 2, median),
	serum=levels(H1N1_mod$serumStrain),
	serumYear=H1N1_mod$serumYear[match(serum, H1N1_mod$serumStrain)]
)

h1n1_virus3 <- data_frame(
	x=apply(ee_h1n1$X[1:2000+4000,1,], 2, median),
	y=apply(ee_h1n1$X[1:2000+4000,2,], 2, median),
	virus=levels(H1N1_mod$virusStrain),
	virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)]
)

h1n1_serum3 <- data_frame(
	x=apply(ee_h1n1$Y[1:2000+4000,1,], 2, median),
	y=apply(ee_h1n1$Y[1:2000+4000,2,], 2, median),
	serum=levels(H1N1_mod$serumStrain),
	serumYear=H1N1_mod$serumYear[match(serum, H1N1_mod$serumStrain)]
)

h1n1_virus4 <- data_frame(
	x=apply(ee_h1n1$X[1:2000+6000,1,], 2, median),
	y=apply(ee_h1n1$X[1:2000+6000,2,], 2, median),
	virus=levels(H1N1_mod$virusStrain),
	virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)]
) 

h1n1_serum4 <- data_frame(
	x=apply(ee_h1n1$Y[1:2000+6000,1,], 2, median),
	y=apply(ee_h1n1$Y[1:2000+6000,2,], 2, median),
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

ggsave("h1n1_map.pdf", gmap, width=8, height=6)

H1N1_rhat <- bayesplot::rhat(H1N1_fit)
hist(H1N1_rhat[xg <- grepl("X", names(H1N1_rhat))])
hist(H1N1_rhat[yg <- grepl("Y", names(H1N1_rhat))])
hist(H1N1_rhat[!(xg | yg)])

distance_data <- lapply(1:8000, function(n) {
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
	mutate(chain=ceiling(as.numeric(iter)/2000)) %>%
	group_by(virusYear, virus, chain) %>%
	summarize(
		mean=mean(dist),
		lwr=quantile(dist, 0.025),
		upr=quantile(dist, 0.975)
	) %>%
	mutate(
		chain=paste0("Chain ", chain)
	)
	
gdist <- ggplot(distance_summary) +
	geom_errorbar(aes(virusYear, ymin=lwr, ymax=upr, group=virus), 
				  width=0, position = position_dodge(width=0.4), col="grey") +
	geom_point(aes(virusYear, mean, group=virus), position = position_dodge(width=0.4)) +
	geom_smooth(aes(virusYear, mean), col="red", se=FALSE) +
	scale_x_continuous("Year") +
	scale_y_continuous("Antigenic distance") +
	facet_wrap(~chain) +
	theme(
		strip.background = element_blank(),
		panel.spacing = grid::unit(0, "cm")
	)

ggsave("h1n1_distance.pdf", gdist, width=8, height=6)

h1n1_virus_effect <- data_frame(
	median=apply(ee_h1n1$re_J, 2, median),
	lwr=apply(ee_h1n1$re_J, 2, quantile, 0.025),
	upr=apply(ee_h1n1$re_J, 2, quantile, 0.975),
	virus=levels(H1N1_mod$virusStrain),
	year=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)],
	key="A/H1N1"
)

h1n1_serum_effect <- data_frame(
	median=apply(ee_h1n1$re_A, 2, median),
	lwr=apply(ee_h1n1$re_A, 2, quantile, 0.025),
	upr=apply(ee_h1n1$re_A, 2, quantile, 0.975),
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

ggsave("h1n1_effect.pdf", gtot, width=6, height=4)

drift_data <- lapply(1:8000, function(n) {
	h1n1_virus1 <- data_frame(
		x=ee_h1n1$X[n,1,],
		y=ee_h1n1$X[n,2,],
		virus=levels(H1N1_mod$virusStrain),
		virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)]
	)
	h1n1_mean_loc <- h1n1_virus1 %>%
		group_by(virusYear) %>%
		summarize(
			x=mean(x), y=mean(y)
		) %>%
		arrange(virusYear) %>%
		filter(virusYear > 2000)
	
	drift <- data.frame(
		year=2003:2009,
		dist=sqrt(diff(h1n1_mean_loc$x)^2 + diff(h1n1_mean_loc$y)^2)
	)
	
	drift
}) %>%
	bind_rows(.id="sim")

drift_data2 <- drift_data  %>%
	mutate(chain=ceiling(as.numeric(sim)/2000)) %>%
	group_by(year, chain) %>%
	summarize(
		mean=mean(dist),
		lwr=quantile(dist, 0.025),
		upr=quantile(dist, 0.975)
	)

avidity_drift <- lapply(1:8000, function(n) {
	dd <- data_frame(
		effect=ee_h1n1$re_J[n,],
		virus=levels(H1N1_mod$virusStrain),
		year=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)],
		key="A/H1N1"
	)  %>%
		filter(year > 2000) %>%
		group_by(year) %>%
		summarize(
			mean=mean(effect)
		)
	
	data.frame(
		year=2003:2009,
		dist2=abs(diff(dd$mean))
	)
})  %>%
	bind_rows(.id="sim")

avidity_drift2 <- avidity_drift %>%
	mutate(chain=ceiling(as.numeric(sim)/2000)) %>%
	group_by(year, chain) %>%
	summarize(
		mean=mean(dist2),
		lwr=quantile(dist2, 0.025),
		upr=quantile(dist2, 0.975)
	)

merge(drift_data, avidity_drift) %>%
	group_by(sim) %>%
	summarize(cor=cor(dist, dist2)) %>%
	summarize(
		mean=mean(cor),
		lwr=quantile(cor, 0.025),
		upr=quantile(cor, 0.975)
	)

gdrift <- ggplot(drift_data2) +
	geom_point(aes(year-0.1, mean, col="Antigenic")) +
	geom_line(aes(year-0.1, mean, col="Antigenic")) +
	geom_errorbar(aes(year-0.1, ymin=lwr, ymax=upr, col="Antigenic"), width=0.1) +
	geom_point(data=avidity_drift2, aes(year, mean*3, col="Non-antigenic")) +
	geom_line(data=avidity_drift2, aes(year, mean*3, col="Non-antigenic")) +
	geom_errorbar(data=avidity_drift2, aes(year, ymin=lwr*3, ymax=upr*3, col="Non-antigenic"), width=0.1) +
	scale_y_continuous("Antigenic drift", 
					   sec.axis = sec_axis(~./3, name = "Non-antigenic drift in virus effect")) +
	scale_x_continuous("Year") +
	facet_wrap(~chain) +
	theme(
		legend.title = element_blank(),
		legend.position = "top"
	)

ggsave("h1n1_drift.pdf", gdrift, width=6, height=4)
