data {
    int E; // number of edges
    int D; // target dimension
    vector<lower=0>[E] HH;
    int<lower=0, upper=1> censoring[E]; // censoring of distances
    /* 0: uncensored, 
     * 1: left-censored
     */
     int nstrain;
     
     int<lower=1, upper=nstrain> alevel[E];
     int<lower=1, upper=nstrain> vlevel[E];
}

parameters {
    /* the first point is fixed at the origin (0)
     * the next D points are stored in a lower triangular matrix (X1)
     * the remaining points are given by a N-D-1 x D matrix (X2)
     */
    cholesky_factor_cov[D] mu1;
    matrix[D, nstrain-D-1] mu2;
    matrix[D, nstrain] X;
    matrix[D, nstrain] Y;
    real<lower=0> sigma; // residual sd
    
    vector[nstrain] JJ; // unscaled virus level effect
    vector[nstrain] AA; // unscaled serum level effect
	
	real<lower=0> sigma_J;
	real<lower=0> sigma_A;
	real intercept;
}

transformed parameters {
    // put all coordinates in a single matrix X
    matrix[D, nstrain] mu = append_col(rep_vector(0, D), append_col(mu1', mu2)); // [0, mu1', mu2]
    
	vector[nstrain] re_J = JJ * sigma_J;
	vector[nstrain] re_A = AA * sigma_A;
	
}

model {
	sigma_J ~ gamma(5, 5);
	sigma_A ~ gamma(5, 5);
	intercept ~ normal(10, 1);
	
    JJ ~ normal(0, 1);
    AA ~ normal(0, 1);
	
	sigma ~ cauchy(0, 10);
	
    for (i in 1:D) {
    	mu1[i,] ~ normal(0, 100);
    	mu2[i,] ~ normal(0, 100);
    	X[i,] ~ normal(mu[i,], 0.1);
    	Y[i,] ~ normal(mu[i,], 0.1);
    }
    
    // calculate the current distances (between the Xs)
    for (e in 1:E) {
        real dist;
        real HHhat;
        
        dist =  distance(col(X, vlevel[e]), col(Y, alevel[e]));
		HHhat = intercept - dist + re_J[vlevel[e]] + re_A[alevel[e]];
        
		// likelihood
        if ( censoring[e] == 0 ) { // uncensored
        	target += log(normal_cdf(HH[e]+1, HHhat, sigma) - normal_cdf(HH[e], HHhat, sigma));
		} else if ( censoring[e] == 1 ) { // left-censored
            target += normal_lcdf(HH[e] | HHhat, sigma);
        }
    }
}
