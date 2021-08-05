// useful functions to use in R

functions {
  vector eta_guess(int[] y,
                   vector ye,
                   vector[] x,
                   real alpha,
                   real rho) {
    int n_obs = size(y);
    vector[n_obs] eta_mle;
    vector[n_obs] theta_mle;

    // hacky way to deal with case where 
    for (i in 1:n_obs) {
      real y_log[2] = {log(y[i]), 0};
      theta_mle[i] = max(y_log) - log(ye[i]);
    }

    print(theta_mle);

    {
      real delta = 1e-8;
      matrix[n_obs, n_obs] L_Sigma;
      matrix[n_obs, n_obs] Sigma;
      Sigma = cov_exp_quad(x, alpha, rho);
      for (n in 1:n_obs) Sigma[n, n] = Sigma[n,n] + delta;
      L_Sigma = cholesky_decompose(Sigma);
      eta_mle = mdivide_left(L_Sigma, theta_mle);
    }
    return eta_mle;
  }
}
