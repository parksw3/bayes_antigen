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
     
     int<lower=1, upper=nvirus> vlevel[E];
     int<lower=1, upper=nserum> slevel[E];
     
}

parameters {
    /* the first point is fixed at the origin (0)
     * the next D points are stored in a lower triangular matrix (X1)
     * the remaining points are given by a N-D-1 x D matrix (X2)
     */
    cholesky_factor_cov[D] X1;
    matrix[D, nvirus-D-1] X2;
    matrix[D, nserum] Y;
    real<lower=0> sigma; // residual sd
    
    vector[nvirus] JJ; // unscaled virus level effect
    
	vector[nserum] AA; // unscaled serum level effect
	
	real<lower=0> sigma_J;
	real<lower=0> sigma_A;
	real intercept;
}

transformed parameters {
    // put all coordinates in a single matrix X
    matrix[D, nvirus] X = append_col(rep_vector(0, D), append_col(X1', X2)); // [0, X1', X2]
    
	vector[nvirus] re_J = JJ * sigma_J;
	vector[nserum] re_A = AA * sigma_A;
	
}

model {
	JJ ~ normal(0, 1);
    AA ~ normal(0, 1);
	
	sigma ~ cauchy(0, 10);
	sigma_J ~ gamma(5, 5);
	sigma_A ~ gamma(5, 5);
	intercept ~ normal(10, 1);
	
    for (i in 1:D) {
    	X1[i,] ~ normal(0, 100);
    	X2[i,] ~ normal(0, 100);
    	Y[i,] ~ normal(0, 100);
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
