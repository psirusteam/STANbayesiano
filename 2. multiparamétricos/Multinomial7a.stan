//this is a multinomial logit Normal Model with correlated random effects
//It reparametrizse the model as in Koster and McElreath (2017), McElreath (2020),
//as well as the Stan User's guide version 2.24, Section 1.13
//Written by Carolina Franco 

///////////////////////////////////////////////////


//Custom multinomial distribution by Joseph Kang 3/8/2021
//This function is to avoid having to round data before entering it in

//Custom multinomial distribution by Joseph Kang 3/8/2021
functions{
  real mtnomial_lpdf(real[] y, vector theta){//note that y is real
    real lprob=0;//hold logprog
    real n = y[1]+y[2]+y[3];
    for(i in 1:3){//3 categories
       lprob = lprob + y[i]*log(theta[i]); //
    }//for
    lprob = lprob + 
    lgamma(n+1)-(lgamma(y[1]+1)+
    lgamma(y[2]+1)+lgamma(y[3]+1));//coefficients were added.
    return lprob;
  } //mtnomial_lpdf
  

}//functions

data{
    int K;
    int D;     //N_obs + N_mis 
    int P;         // number of possible categories
    real y[D,P];   // MAKE SURE NO NA'S FOR MISSING DATA: STAN DOES NOT ALLOW 
    matrix [D, K] X; // X[D] row_vector of size D;

//NOTE:  X, X, should include INTERCEPT

}//data



parameters{
    matrix[P-1,D]                 z;//matrix of standardized random effects
    vector<lower=0>[P-1]      sigma;//standard deviation of random effects
    cholesky_factor_corr[P-1] L_Rho;//Cholesky factor of correlation matrix of random effects
    vector[K] beta1;
    vector[K] beta2; //regression coefficients
}//parameters

transformed parameters{
    matrix[D, P-1]      u;     // matrix of scaled random effects
    matrix[D, P-1]      p;     // conditional probabilities (see below)
    matrix[D, P]        theta; // prob of disjoint categories (see below)
         
	u     = (diag_pre_multiply(sigma,L_Rho)*z)' ;// note transpose in this transformation
    	
	for(i in 1:D){
	
	p[i, 1]=inv_logit(u[i, 1]+X[i, ]*beta1); 	// p1= any hearing loss
	p[i, 2]=inv_logit(u[i, 2]+X[i, ]*beta2); 	// p2= severe given any hearing loss

      
        theta[i, 1]=1-p[i, 1];                  	// theta1=normal hearing prob 
        theta[i, 2]=p[i, 1]*(1-p[i, 2]);         	// theta2=mild but not severe hearing loss prob
        theta[i, 3]=p[i, 1]*p[i, 2];    		// theta3=severe hearing loss prob	
    }//i

}//transformed parameters

model{
    to_vector(z)     ~ normal(0,1);
    to_vector(beta1) ~ normal(0,100); 
    to_vector(beta2) ~ normal(0,100); 
    sigma ~ cauchy(0,2.5);//exponential(1);
    L_Rho ~ lkj_corr_cholesky(1); 
    
    for (i in 1:D){ 
	y[i, ] ~ multinomial(theta[i, ]); 
    }//i
}//model

generated quantities{
    matrix[P-1,P-1] Omega; //The correlation matrix   
    matrix[P-1,P-1] Sigma; //The covariance matrix 
	vector[D] log_lik;  //for computing lpmf and 

    Omega = multiply_lower_tri_self_transpose(L_Rho);
    Sigma = quad_form_diag(Omega,sigma);
	for (i in 1:D)
	log_lik[i] = mtnomial_lpdf(y[i,] |theta[ ,i]);//loglikelihood
}//generatedquantities





