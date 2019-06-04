library(dplyr)
library(bayesplot)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

load("../analysis/h3n2_analysis.rda")
load("../analysis/h1n1_analysis.rda")
load("../analysis/victoria_analysis.rda")
load("../analysis/yamagata_analysis.rda")

FluNetdata <- read.csv("../data/FluNetInteractiveReport_NA_1995_2011.csv", skip=2)

flusumm <- FluNetdata %>%
	group_by(Year, Country) %>%
	summarize(
		AH3=sum(AH3, na.rm=TRUE)/sum(SPEC_PROCESSED_NB, na.rm=TRUE),
		AH1=(sum(AH1, na.rm=TRUE) + sum(AH1N12009, na.rm=TRUE))/sum(SPEC_PROCESSED_NB, na.rm=TRUE)
	) %>%
	summarize(
		AH3=mean(AH3),
		AH1=mean(AH1)
	) %>%
	filter(Year <= 2011)

ee_h3n2 <- rstan::extract(H3N2_fit)
ee_h1n1 <- rstan::extract(H1N1_fit)
ee_vic <- rstan::extract(Vic_fit)
ee_yam <- rstan::extract(Yam_fit)

h3n2_effect <- data_frame(
	median=apply(apply(ee_h3n2$JJ, 2, "*", ee_h3n2$sigma_J), 2, median),
	lwr=apply(apply(ee_h3n2$JJ, 2, "*", ee_h3n2$sigma_J), 2, quantile, 0.025),
	upr=apply(apply(ee_h3n2$JJ, 2, "*", ee_h3n2$sigma_J), 2, quantile, 0.975),
	virus=levels(H3N2_mod$virusStrain),
	virusYear=H3N2_mod$virusYear[match(virus, H3N2_mod$virusStrain)],
	key="A/H3N2"
)

h1n1_effect <- data_frame(
	median=apply(apply(ee_h1n1$JJ, 2, "*", ee_h1n1$sigma_J), 2, median),
	lwr=apply(apply(ee_h1n1$JJ, 2, "*", ee_h1n1$sigma_J), 2, quantile, 0.025),
	upr=apply(apply(ee_h1n1$JJ, 2, "*", ee_h1n1$sigma_J), 2, quantile, 0.975),
	virus=levels(H1N1_mod$virusStrain),
	virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)],
	key="A/H1N1"
)

vic_effect <- data_frame(
	median=apply(apply(ee_vic$JJ, 2, "*", ee_vic$sigma_J), 2, median),
	lwr=apply(apply(ee_vic$JJ, 2, "*", ee_vic$sigma_J), 2, quantile, 0.025),
	upr=apply(apply(ee_vic$JJ, 2, "*", ee_vic$sigma_J), 2, quantile, 0.975),
	virus=levels(Vic_mod$virusStrain),
	virusYear=Vic_mod$virusYear[match(virus, Vic_mod$virusStrain)],
	key="B/Vic"
)

yam_effect <- data_frame(
	median=apply(apply(ee_yam$JJ, 2, "*", ee_yam$sigma_J), 2, median),
	lwr=apply(apply(ee_yam$JJ, 2, "*", ee_yam$sigma_J), 2, quantile, 0.025),
	upr=apply(apply(ee_yam$JJ, 2, "*", ee_yam$sigma_J), 2, quantile, 0.975),
	virus=levels(Yam_mod$virusStrain),
	virusYear=Yam_mod$virusYear[match(virus, Yam_mod$virusStrain)],
	key="B/Yam"
)

JJ_effect <- rbind(h3n2_effect, h1n1_effect, vic_effect, yam_effect)

gavid <- ggplot(JJ_effect) +
	geom_hline(yintercept=0, lty=2) +
	geom_boxplot(aes(virusYear, median, group=interaction(virusYear), fill=key, col=key), alpha=0.1) +
	geom_smooth(aes(virusYear, median, col=key), span=0.2, se=FALSE, lwd=1.5, n=300) +
	scale_y_continuous("Random effects of virus strain") +
	scale_x_continuous("Year") +
	facet_grid(key~., scale="free_y") +
	theme(
		panel.grid = element_blank(),
		panel.spacing = grid::unit(0, "cm"),
		legend.position = "none",
		strip.background = element_blank()
	)

ggsave("avidity_summary.pdf", gavid, width=6, height=4)

h3n2_cor <- h3n2_effect %>%
	filter(virusYear %in% flusumm$Year) %>%
	mutate(Year=virusYear) %>%
	merge(flusumm)

h3n2_cor_eff <- apply(ee_h3n2$JJ, 1, function(x) {
	dd <- data_frame(
		re=x,
		virus=levels(H3N2_mod$virusStrain),
		virusYear=H3N2_mod$virusYear[match(virus, H3N2_mod$virusStrain)]
	)  %>%
		filter(virusYear %in% flusumm$Year) %>%
		mutate(Year=virusYear) %>%
		merge(flusumm)
	
	cor(dd$re, dd$AH3)
})

h1n1_cor <- h1n1_effect %>%
	filter(virusYear %in% flusumm$Year) %>%
	mutate(Year=virusYear) %>%
	merge(flusumm)

h1n1_cor_eff <- apply(ee_h1n1$JJ, 1, function(x) {
	dd <- data_frame(
		re=x,
		virus=levels(H1N1_mod$virusStrain),
		virusYear=H1N1_mod$virusYear[match(virus, H1N1_mod$virusStrain)]
	)  %>%
		filter(virusYear %in% flusumm$Year) %>%
		mutate(Year=virusYear) %>%
		merge(flusumm)
	
	cor(dd$re, dd$AH1)
})

g1 <- ggplot(h3n2_cor) +
	geom_point(aes(AH3, median)) +
	geom_smooth(aes(AH3, median), method="lm", col="black", lty=2) +
	annotate("text", x=0.035, y=3.7, 
			 label=paste0("Correlation coefficient: ",
			 			 round(median(h3n2_cor_eff), 2), 
			 			 " (",
			 			 paste(round(quantile(h3n2_cor_eff, c(0.025, 0.975)), 2), collapse=", ")
			 			 , ")")) +
	scale_x_continuous("Proportion of positive A3N2 cases (1997-2011)") +
	scale_y_continuous("Random effects of virus strain")

g2 <- ggplot(h1n1_cor) +
	geom_point(aes(AH1, median)) +
	geom_smooth(aes(AH1, median), method="lm", col="black", lty=2) +
	annotate("text", x=0.092, y=2, 
			 label=paste0("Correlation coefficient: ",
			 			 round(median(h1n1_cor_eff), 2), 
			 			 " (",
			 			 paste(round(quantile(h1n1_cor_eff, c(0.025, 0.975)), 2), collapse=", ")
			 			 , ")")) +
	scale_x_continuous("Proportion of positive A1N1 cases (1999-2009)") +
	scale_y_continuous("Random effects of virus strain")

gbind <- arrangeGrob(g1, g2, nrow=1)

ggsave("correlation.pdf", gbind, width=10, height=4)
