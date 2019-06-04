library(dplyr)
library(rstan)

H1N1_data <- read.table("../data/H1N1_HI_data.tsv", sep = '\t', header = TRUE)

H1N1_mod <- H1N1_data %>%
	arrange(virusYear, serumYear) %>%
	mutate(
		HH=log2(as.numeric(gsub("<", "", titer))),
		censor=grepl("<", titer)
	)

virus_effect <- H1N1_mod %>%
	group_by(virusStrain) %>%
	summarize(max=max(HH)) %>%
	summarize(
		mean=mean(max),
		var=var(max)
	)

serum_effect <- H1N1_mod %>%
	group_by(serumStrain) %>%
	summarize(max=max(HH)) %>%
	summarize(
		mean=mean(max),
		var=var(max)
	)

H1N1_standata <- list(
	E=nrow(H1N1_mod),
	D=2,
	nvirus=length(unique(H1N1_mod$virusStrain)),
	nserum=length(unique(H1N1_mod$serumStrain)),
	vlevel=as.numeric(factor(H1N1_mod$virusStrain)),
	slevel=as.numeric(factor(H1N1_mod$serumStrain)),
	HH=H1N1_mod$HH,
	censoring=as.numeric(H1N1_mod$censor),
	mu_J=virus_effect$mean,
	sigma_J=virus_effect$var,
	mu_A=serum_effect$mean,
	sigma_A=serum_effect$var
)

rt <- stanc(file="../stancode/antigen_model_empirical_bayes.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

H1N1_fit <- sampling(sm, data=H1N1_standata, chain=4, iter=4000, seed=101,
					 control=list(max_treedepth = 15))

save("H1N1_fit", "H1N1_mod", "H1N1_standata", file="h1n1_analysis.rda")
