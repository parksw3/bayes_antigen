library(rstan)
library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(gridExtra)

load("../figure/hamming_map.rda")

load("../analysis/h1n1_hamming_informed_analysis.rda")

ee <- rstan::extract(H1N1_fit)

ee_rhat <- bayesplot::rhat(H1N1_fit)

## converged!!!
hist(ee_rhat)
max(ee_rhat)

medianpos <- data.frame(
	x=apply(ee$X[,1,], 2, median),
	y=apply(ee$X[,2,], 2, median),
	name=combined_map[H1N1_standata$vlevel_all,]$name,
	year=combined_map[H1N1_standata$vlevel_all,]$year
)

twomap <- merge(medianpos, combined_map)

g1 <- ggplot(combined_map) +
	geom_point(aes(V1, V2, col=year)) +
	scale_color_gradient(low="red", high="blue") +
	xlab("Dimension 1") +
	ylab("Dimension 2")
	
g2 <- ggplot(medianpos) +
	geom_point(aes(x, y, col=year)) +
	scale_color_gradient(low="red", high="blue") +
	xlab("Dimension 1") +
	ylab("Dimension 2")
	
g3 <- ggplot(twomap) +
	geom_point(aes(V1, V2, col="sequence"), size=2) +
	geom_point(aes(x, y, col="titer"), size=2) +
	geom_segment(aes(V1, V2, xend=x, yend=y))  +
	xlab("Dimension 1") +
	ylab("Dimension 2")

virus_effect <- data.frame(
	median=apply(apply(ee$JJ, 2, "*", ee$sigma_J), 2, median),
	lwr=apply(apply(ee$JJ, 2, "*", ee$sigma_J), 2, quantile, 0.025),
	upr=apply(apply(ee$JJ, 2, "*", ee$sigma_J), 2, quantile, 0.975),
	year=combined_map[H1N1_standata$vlevel_all,]$year,
	name=combined_map[H1N1_standata$vlevel_all,]$name
)

g4 <- ggplot(virus_effect) +
	geom_point(aes(year, median, group=name), position=position_dodge(0.5)) +
	geom_errorbar(aes(year, ymin=lwr, ymax=upr, group=name), width=0, position=position_dodge(0.5)) +
	geom_smooth(aes(year, median), span=0.4, col="red", se=FALSE)

serum_effect <- data.frame(
	median=apply(apply(ee$AA, 2, "*", ee$sigma_A), 2, median),
	lwr=apply(apply(ee$AA, 2, "*", ee$sigma_A), 2, quantile, 0.025),
	upr=apply(apply(ee$AA, 2, "*", ee$sigma_A), 2, quantile, 0.975),
	year=combined_map[H1N1_standata$slevel_all,]$year,
	name=combined_map[H1N1_standata$slevel_all,]$name
)

g5 <- g4 %+% serum_effect

gcombine1 <- arrangeGrob(
	g1 + ggtitle("Sequence only"),
	g2 + ggtitle("Sequence + HI titer"),
	g3 + ggtitle("Sequence + HI titer (direction)"),
	nrow=1
)

gcombine2 <- arrangeGrob(
	g4 + ylab("Virus effect"),
	g5 + ylab("Serum effect"),
	nrow=2
)

ggsave("h1n1_hamming_informed_map.pdf", gcombine1, width=12, height=4)
ggsave("h1n1_hamming_informed_effect.pdf", gcombine2, width=6, height=4)
