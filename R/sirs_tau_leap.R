SIRS_tau_leap <- function(population,   # vector length n_locs
                          sir,          # matrix [n_locs, 3]: S, I, R
                          beta,         # matrix [n_locs, 1]
                          t_steps = 7,
                          tau     = 1,
                          gamma   = 1/7,
                          omega   = 1/180) {
  
  if (is.null(dim(sir))) sir <- matrix(sir, nrow = 1)
  n_locs <- nrow(sir)
  
  out <- vector("list", t_steps + 1)
  out[[1]] <- sir
  
  state <- sir
  for (step in seq_len(t_steps)) {
    S <- state[, 1]; I <- state[, 2]; R <- state[, 3]
    beta_vec <- as.numeric(beta)
    
    lambda_SI <- pmax(beta_vec * S * I / population, 0) * tau
    lambda_IR <- pmax(gamma * I, 0) * tau
    lambda_RS <- pmax(omega * R, 0) * tau
    
    dSI <- rpois(n_locs, lambda_SI)
    dIR <- rpois(n_locs, lambda_IR)
    dRS <- rpois(n_locs, lambda_RS)
    
    # Cap so no compartment goes negative
    dSI <- pmin(dSI, S)
    dIR <- pmin(dIR, I)
    dRS <- pmin(dRS, R)
    
    new_state <- matrix(0, nrow = n_locs, ncol = 3)
    new_state[, 1] <- S - dSI + dRS
    new_state[, 2] <- I + dSI - dIR
    new_state[, 3] <- R + dIR - dRS
    
    state <- new_state
    out[[step + 1]] <- state
  }
  
  result <- array(0, dim = c(n_locs, 3, t_steps + 1))
  for (s in seq_along(out)) {
    result[, , s] <- out[[s]]
  }
  result
}