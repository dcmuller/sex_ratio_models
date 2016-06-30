## poisson model for rates
data {
  int<lower=0> N;
  int<lower=1> n_coef;
  int<lower=0> y[N];
  vector[N] log_offset;
  int<lower=1> n_reg;
  int<lower=1, upper=n_reg> reg[N];
  int<lower=1> n_age;
  int<lower=1, upper=n_age> age[N];
  int<lower=1> n_year;
  int<lower=1, upper=n_year> year[N];
  row_vector[n_coef] x[N];
}
parameters {
  real log_lambda; //average log incidence
  real<lower=0> sigma_log_lambda;

  vector[n_age] age_intercept; // log baseline rate by age group
  real<lower=0> sigma_age;

  vector[n_year] year_intercept; // year-specific deviation from log baseline rate
  real<lower=0> sigma_year; // sd of cohort log baseline rates

  vector[n_reg] reg_intercept; // year-specific deviation from log baseline rate
  real<lower=0> sigma_reg; // sd of cohort log baseline rates

  real<lower=0> theta; // alpha and beta parameters of gamma distribution for overdispersion
  real<lower=0> h; // beta parameter of gamma prior on theta
  real<lower=0> od[N]; // overdispersion parameters

  vector[n_coef] mu_b; // averave regression coefs
  vector[n_coef] b_age_raw[n_age]; // array of vectors for regression coefs
  vector[n_coef] b_year_raw[n_year]; // array of vectors for regression coefs
  vector[n_coef] b_reg_raw[n_reg]; // array of vectors for regression coefs
  real<lower=0> sigma_b_age[n_coef]; // sds of regression coefs
  real<lower=0> sigma_b_year[n_coef]; // sds of regression coefs
  real<lower=0> sigma_b_reg[n_coef]; // sds of regression coefs
}
transformed parameters {
  vector[n_coef] b_age[n_age]; // array of vectors for regression coefs
  vector[n_coef] b_year[n_year]; // array of vectors for regression coefs
  vector[n_coef] b_reg[n_reg]; // array of vectors for regression coefs
  for (a in 1:n_age) {
    for (j in 1:n_coef)
      b_age[a][j] <- b_age_raw[a][j] * sigma_b_age[j];
  }
  for (i in 1:n_year) {
    for (j in 1:n_coef)
      b_year[i][j] <- b_year_raw[i][j] * sigma_b_year[j];
  }
  for (r in 1:n_reg) {
    for (j in 1:n_coef)
      b_reg[r][j] <- b_reg_raw[r][j] * sigma_b_reg[j];
  }
}
model {
  vector[N] lambda_hat_log;
  // priors
  sigma_log_lambda ~ normal(0,2);
  log_lambda ~ normal(-9, sigma_log_lambda);

  sigma_age ~ normal(0,1); // sd of average rate (by age)
  age_intercept ~ normal(0, sigma_age); // average rate by year

  sigma_year ~ normal(0, 0.2);         // weakly informative for between registry sd
  year_intercept ~ normal(0, sigma_year);
  sigma_reg ~ normal(0, 0.2);         // weakly informative for between registry sd
  reg_intercept ~ normal(0, sigma_reg);

  h ~ gamma(0.01, 0.01); // hyperprior for beta parameter of gamma distributed overdispersion
  theta ~ gamma(0.01, h); // alpha and beta parameter of overdispersion

  mu_b ~ normal(0, 2);                   // weakly informative prior for log rate ratios
  sigma_b_age ~ normal(0, 0.2);           // weakly informative for sd of regression coefs
  sigma_b_year ~ normal(0, 0.2);           // weakly informative for sd of regression coefs
  sigma_b_reg ~ normal(0, 0.2);           // weakly informative for sd of regression coefs

  // likelihood
  od ~ gamma(theta, theta); // overdispersion
  for (a in 1:n_age)
    b_age_raw[a] ~ normal(0, 1); // normal model for regression coefs
  for (i in 1:n_year)
    b_year_raw[i] ~ normal(0, 1); // normal model for regression coefs
  for (r in 1:n_reg)
    b_reg_raw[r] ~ normal(0, 1); // normal model for regression coefs
  // predicted rate
  for (n in 1:N)
    lambda_hat_log[n] <- log_lambda + age_intercept[age[n]] +
      year_intercept[year[n]] + reg_intercept[reg[n]] +
      x[n] * (mu_b + b_age[age[n]] + b_year[year[n]] + b_reg[reg[n]]) +
      log_offset[n] + log(od[n]);

  increment_log_prob(poisson_log_log(y, lambda_hat_log));
}
generated quantities {
  vector[n_coef] log_rr_age[n_age];
  vector[n_coef] log_rr_year[n_year];
  vector[n_coef] log_rr_reg[n_reg];
  vector[n_coef] log_sigma_ratio_age;
  vector[n_coef] log_sigma_ratio_year;
  vector[n_coef] log_sigma_ratio_reg;

  for (i in 1:n_age)
    log_rr_age[i] <- mu_b + b_age[i];
  for (i in 1:n_year)
    log_rr_year[i] <- mu_b + b_year[i];
  for (i in 1:n_reg)
    log_rr_reg[i] <- mu_b + b_reg[i];

  log_sigma_ratio_age <- log(to_vector(sigma_b_age)) - log(sigma_age);
  log_sigma_ratio_year <- log(to_vector(sigma_b_year)) - log(sigma_year);
  log_sigma_ratio_reg <- log(to_vector(sigma_b_reg)) - log(sigma_reg);

#  real log_baseline_rate[n_age, n_year, n_reg];
#  vector[n_coef] log_rr[n_age, n_year, n_reg];
#  for (a in 1:n_age) {
#    for (i in 1:n_year) {
#      for (r in 1:n_reg) {
#        log_baseline_rate[a, i, r] <- log_lambda0[a] + year_intercept[i] + reg_intercept[r];
#        log_rr[a, i, r] <- mu_b + b_age[a] + b_year[i] + b_reg[r];
#      }
#    }
#  }
}
