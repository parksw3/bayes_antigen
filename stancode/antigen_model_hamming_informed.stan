data {
    int E; // number of edges
    int D; // target dimension
    vector<lower=0>[E] HH;
    int<lower=0, upper=1> censoring[E]; // censoring of distances
    /* 0: uncensored, 
     * 1: left-censored
     */
     int nvirus;
     int nserum;
     int nall; // number of all strains
     
     int<lower=1, upper=nvirus> vlevel[E];
     int<lower=1, upper=nserum> slevel[E];
     
     int<lower=1, upper=nall> vlevel_all[nvirus];
     int<lower=1, upper=nall> slevel_all[nserum];
     
     matrix[D, nall] prior_mean;
}

parameters {
    matrix[D, nvirus] X;
    matrix[D, nserum] Y;
    real<lower=0> sigma; // residual sd
    
    vector[nvirus] JJ; // unscaled virus level effect
    vector[nserum] AA; // unscaled serum level effect
	
	real<lower=0> sigma_J;
	real<lower=0> sigma_A;
	real intercept;
}

transformed parameters {
    vector[nvirus] re_J = JJ * sigma_J;
	vector[nserum] re_A = AA * sigma_A;
}

model {
	sigma_J ~ gamma(5, 5);
	sigma_A ~ gamma(5, 5);
	intercept ~ normal(0, 5);
	
    JJ ~ normal(0, 1);
    AA ~ normal(0, 1);
	
	sigma ~ cauchy(0, 10);
	
	for (i in 1:nvirus) {
		for (j in 1:D) {
			X[j,i] ~ normal(prior_mean[j,vlevel_all[i]], sqrt(0.5));
		}
	}
	
	for (i in 1:nserum) {
		for (j in 1:D) {
			Y[j,i] ~ normal(prior_mean[j,slevel_all[i]], sqrt(0.5));
		}
	}
	
	// calculate the current distances (between the Xs)
    for (e in 1:E) {
        real dist;
        real HHhat;
        
        dist =  distance(col(X, vlevel[e]), col(Y, slevel[e]));
		HHhat = intercept - dist + re_J[vlevel[e]] + re_A[slevel[e]];
        
		// likelihood
        if ( censoring[e] == 0 ) { // uncensored
        	target += log(normal_cdf(HH[e]+1, HHhat, sigma) - normal_cdf(HH[e], HHhat, sigma));
		} else if ( censoring[e] == 1 ) { // left-censored
            target += normal_lcdf(HH[e] | HHhat, sigma);
        }
    }
}
