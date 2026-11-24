# ---------------------------------------------------------------------------
# One-step forecast: SIRS tau-leap + SEASONAL beta OU process
# ---------------------------------------------------------------------------
# sir_particles:  array [n_locs, 3, N]
# beta_particles: matrix [n_locs, N]
# current_date:   Date object; used to compute seasonal mean of log(beta)
#
# Seasonal OU on log(beta):
#   log_mu(t) = log(beta0) + a * cos(2*pi*(doy - doy_peak)/365)
#   d log(beta) = theta*(log_mu(t) - log(beta)) dt + sigma dW
# ---------------------------------------------------------------------------

.doy <- function(date) as.integer(format(date, "%j"))

seasonal_log_mu <- function(date, beta0, a, doy_peak) {
  doy <- .doy(date)
  log(beta0) + a * cos(2 * pi * (doy - doy_peak) / 365)
}

forecast_step <- function(sir_particles,
                          beta_particles,
                          population,
                          current_date,
                          steps_per_obs = 7,
                          theta    = 1/200,
                          beta0    = 0.25,
                          a_season = 0.4, 
                          doy_peak = 25,
                          sigma    = 0.06,
                          gamma = 1/7,
                          omega = 1/180) {
  dims   <- dim(beta_particles)
  n_locs <- dims[1]
  N      <- dims[2]
  
  sir_next <- array(0, dim = dim(sir_particles))
  
  # 1. State forecast: SIRS tau-leap for each particle
  for (p in seq_len(N)) {
    sir_slice <- sir_particles[, , p, drop = FALSE]
    sir_mat   <- matrix(sir_slice, nrow = n_locs, ncol = 3)
    
    beta_p <- matrix(beta_particles[, p], nrow = n_locs, ncol = 1)
    
    traj <- SIRS_tau_leap(population = population,
                          sir        = sir_mat,
                          beta       = beta_p,
                          t_steps    = steps_per_obs,
                          gamma = gamma, omega = omega)
    
    last_t <- dim(traj)[3]
    traj_slice <- traj[, , last_t, drop = FALSE]
    sir_next[, , p] <- matrix(traj_slice, nrow = n_locs, ncol = 3)
  }
  sir_next[sir_next < 0] <- 0
  
  # 2. Beta OU forecast with seasonal forcing
  dt        <- steps_per_obs
  mid_date  <- current_date + floor(dt / 2)
  log_mu_t  <- seasonal_log_mu(mid_date, beta0, a_season, doy_peak)
  
  A <- exp(-theta * dt)
  M <- log_mu_t * (1 - A)
  C <- sigma * sqrt((1 - A^2) / (2 * theta))
  
  noise <- matrix(rnorm(n_locs * N), nrow = n_locs, ncol = N)
  log_beta_next <- A * log(beta_particles) + M + C * noise
  beta_next     <- exp(log_beta_next)
  
  list(sir = sir_next, beta = beta_next)
}