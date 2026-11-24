# ---------------------------------------------------------------------------
# Perturbed-observation EnKF update for SIRS model
# ---------------------------------------------------------------------------
# Joint state per location: [S, I, R, beta] (4 dims)
# Observation: I compartment only (row 2 of the joint state)
# Locations are independent; we run a separate EnKF per location.
# ---------------------------------------------------------------------------

enkf_update_step <- function(sir_init,      # [n_locs, 3, N]
                             beta_init,     # [n_locs, N]
                             observation,   # vector length n_locs (observed I counts)
                             fixed_sigma) { # observation noise std (in count units)
  dims   <- dim(sir_init)
  n_locs <- dims[1]
  N      <- dims[3]
  n_vars <- 4  # S, I, R, beta
  
  sir_post  <- array(0, dim = dim(sir_init))
  beta_post <- matrix(0, nrow = n_locs, ncol = N)
  
  for (loc in seq_len(n_locs)) {
    # Build joint ensemble for this location: [4, N]
    sir_slice <- sir_init[loc, , , drop = FALSE]          # [1, 3, N]
    sir_mat   <- matrix(sir_slice, nrow = 3, ncol = N)    # [3, N]
    ensemble  <- rbind(sir_mat, beta_init[loc, ])         # [4, N]
    
    # Observation operator: picks out I (row 2)
    H_op <- matrix(0, nrow = 1, ncol = n_vars); H_op[1, 2] <- 1
    R_obs <- matrix(fixed_sigma^2, 1, 1)
    
    prior_cov <- cov(t(ensemble))                         # [4, 4]
    HPHt      <- H_op %*% prior_cov %*% t(H_op)           # scalar
    K         <- prior_cov %*% t(H_op) %*% solve(R_obs + HPHt)  # [4, 1]
    
    # Perturbed observations
    eps   <- rnorm(N, 0, fixed_sigma)
    y     <- observation[loc]
    innov <- (y - ensemble[2, ]) + eps                    # length N
    post  <- ensemble + K %*% matrix(innov, nrow = 1)
    
    # Split back out; enforce non-negativity
    sir_post_mat <- post[1:3, , drop = FALSE]
    sir_post_mat[sir_post_mat < 0] <- 0
    sir_post[loc, , ] <- sir_post_mat
    
    beta_vec <- post[4, ]
    beta_vec[beta_vec < 1e-4] <- 1e-4
    beta_post[loc, ] <- beta_vec
  }
  
  list(sir = sir_post, beta = beta_post)
}