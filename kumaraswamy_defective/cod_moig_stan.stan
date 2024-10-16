
  data {
   int<lower=0> N;                        // tamanho amostral
   vector[N] time;                        // tempo de falha observado
   array[N] int<lower=0, upper=1> delta;        // indicador do evento
  }
  
  parameters {
    real alpha;
    real<lower=0> beta;
    real<lower=0> lambda;
  }
  
  model {
    
    alpha~normal(0,10);
    beta~gamma(0.001,0.001);
    lambda~gamma(0.001,0.001);
    
    // Definição manual da função de verossimilhança
    for (i in 1:N) {
      real z1 = (-1 + alpha * time[i]) / sqrt(beta * time[i]);
      real z2 = (-1 - alpha * time[i]) / sqrt(beta * time[i]);
      real phi1 = Phi(z1);
      real phi2 = Phi(z2);
      real exp_term = exp(2 * alpha / beta);

    // Verossimilhança do evento observado e censurado
	
	  target += delta[i] * log(lambda*exp((-1/(2*beta*time[i]))*(1-alpha*time[i])^2)) -
              delta[i] * log(sqrt(2*pi()*beta*time[i]^3) * (lambda + (1-lambda)*(phi1 + exp(2*alpha/beta)*phi2))^2) +
              (1-delta[i])*log(lambda-lambda*(phi1 + exp(2*alpha/beta)*phi2)) -
              (1-delta[i])*log(lambda+(1-lambda)*(phi1 + exp(2*alpha/beta)*phi2));
    }
  }

