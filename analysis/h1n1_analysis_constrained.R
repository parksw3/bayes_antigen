library(dplyr)
library(rstan)

H1N1_data <- read.table("../data/H1N1_HI_data.tsv", sep = '\t', header = TRUE)

H1N1_mod <- H1N1_data %>%
	arrange(virusYear, serumYear) %>%
	mutate(
		HH=log2(as.numeric(gsub("<", "", titer))),
		censor=grepl("<", titer)
	)

allstrain <- unique(c(levels(H1N1_data$virusStrain), levels(H1N1_data$serumStrain)))

H1N1_standata <- list(
	E=nrow(H1N1_mod),
	D=2,
	nstrain=length(allstrain),
	nserum=length(unique(H1N1_mod$serumStrain)),
	vlevel=as.numeric(factor(H1N1_mod$virusStrain, levels=allstrain)),
	alevel=as.numeric(factor(H1N1_mod$serumStrain, levels=allstrain)),
	HH=H1N1_mod$HH,
	censoring=as.numeric(H1N1_mod$censor)
)

rt <- stanc(file="../stancode/antigen_model_constrained.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)

H1N1_fit <- sampling(sm, data=H1N1_standata, chain=1, iter=1000, seed=101,
					 control=list(max_treedepth = 15))

save("H1N1_fit", "H1N1_mod", "H1N1_standata", file="h1n1_analysis_constrained.rda")

oo <- optimizing(sm, data=H1N1_standata, as_vector=FALSE, seed=101,
				 iter=3000)

dist <- sqrt(rowSums((t(oo$par$X)[H1N1_standata$vlevel,] - t(oo$par$Y)[H1N1_standata$alevel,])^2))

plot(oo$par$intercept + oo$par$re_J[H1N1_standata$vlevel] -dist + oo$par$re_A[H1N1_standata$alevel], H1N1_standata$HH)
abline(a=0, b=1, col=2)

plot(t(oo$par$mu))


