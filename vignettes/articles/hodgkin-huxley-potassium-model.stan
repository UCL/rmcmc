// Hodgkin-Huxley voltage clamp model for potassium channel conductance
functions {
  vector alpha_k(vector v, vector k_alpha) {
    return k_alpha[1] .* (v + k_alpha[2])
           ./ (exp((v + k_alpha[2]) ./ k_alpha[3]) - 1);
  }
  vector beta_k(vector v, vector k_beta) {
    return k_beta[1] .* exp(v ./ k_beta[2]);
  }
  vector n_infty(vector v, vector k_alpha, vector k_beta) {
    return alpha_k(v, k_alpha) ./ (alpha_k(v, k_alpha) + beta_k(v, k_beta));
  }
  vector tau_n(vector v, vector k_alpha, vector k_beta) {
    return 1 ./ (alpha_k(v, k_alpha) + beta_k(v, k_beta));
  }
  vector potassium_channel_subunit_activation(
    vector t, vector v, vector k_alpha, vector k_beta
  ) {
    real n_0 = n_infty([0.0]', k_alpha, k_beta)[1]; // Assume starting from equilibrium
    return n_0
           + (n_infty(v, k_alpha, k_beta) - n_0)
             .* (1 - exp(-t ./ tau_n(v, k_alpha, k_beta)));
  }
  vector potassium_conductance(
    vector t, vector v, real g_bar_k, vector k_alpha, vector k_beta, int N
  ) {
    vector[N] n = potassium_channel_subunit_activation(t, v, k_alpha, k_beta);
    return g_bar_k .* n .^ 4;
  }
}
data {
  int<lower=1> N;
  vector[N] times;
  vector[N] depolarizations;
  vector[N] conductances;
}
parameters {
  vector<lower=0>[3] k_alpha;
  vector<lower=0>[2] k_beta;
  real<lower=0> g_bar_k;
  real<lower=0> sigma;
}
model {
  //priors
  k_alpha[1] ~ lognormal(-3, 1);
  k_alpha[2] ~ lognormal(2, 1);
  k_alpha[3] ~ lognormal(2, 1);
  k_beta[1] ~ lognormal(-3, 1);
  k_beta[2] ~ lognormal(2, 1);
  g_bar_k ~ lognormal(2, 1);
  sigma ~ lognormal(0, 1);
  // Solve model to simulate conductances
  vector[N] simulated_conductances = potassium_conductance(
    times, depolarizations, g_bar_k, k_alpha, k_beta, N
  );
  // likelihood
  conductances ~ normal(simulated_conductances, sigma);
}
generated quantities {
  vector[N] simulated_conductances = potassium_conductance(
    times, depolarizations, g_bar_k, k_alpha, k_beta, N
  );
}
