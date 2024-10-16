
  data {
   int<lower=0> N;                        // tamanho amostral
   vector[N] time;                        // tempo de falha observado
   array[N] int<lower=0, upper=1> delta;        // indicador do evento
   real beta;
   real lambda;
  }
  
  parameters {
    real alpha;
    //real<lower=0> beta;
    //real<lower=0> lambda;
  }
  
  model {
    // Prior não informativa para lambda
    
    alpha~normal(0,10);
    //beta~gamma(0.001,0.001);
    //lambda~cauchy(0,5);
    
    // Definição manual da função de verossimilhança
    for (i in 1:N) {
      real z1 = (-1 + alpha * time[i]) / sqrt(beta * time[i]);
      real z2 = (-1 - alpha * time[i]) / sqrt(beta * time[i]);
      real z3 = (1 + alpha * time[i]) / sqrt(beta * time[i]);
      real phi1 = Phi(z1);
      real phi2 = Phi(z2);
      real phi3 = Phi(z3);
      real exp_term = exp(2 * alpha / beta);

    // Verossimilhança do evento observado e censurado
	
	  target += delta[i]*log((lambda*exp(-(alpha*time[i]-1)^2/2*beta*time[i]))/(sqrt(beta*time[i]^3)*((lambda-1)*phi1 + (lambda-1)*exp_term*phi2 - lambda)^2)) +
        		 (1-delta[i])*log((lambda*(phi1 + exp_term*phi2 - 1))/(-(lambda-1)*phi1+(lambda-1)*exp_term*(phi3+1)-1));
    }
  }

