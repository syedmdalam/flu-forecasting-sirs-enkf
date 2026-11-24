# ---------------------------------------------------------------------------
# SIHR tau-leap stochastic simulator (independent locations)
# ---------------------------------------------------------------------------
# State per location: S, I, H, R (counts)
# Transitions (rates):
#   S -> I : beta * S * I / N
#   I -> H : gamma * alpha * I
#   I -> R : gamma * (1 - alpha) * I
#   H -> R : eta * H
# Each transition count per tau step is drawn from Poisson(rate * tau),
# then bounded so no compartment goes negative.
# ---------------------------------------------------------------------------

SIHR_tau_leap <- function(population,   # vector length n_locs
                          sihr,         # matrix [n_locs, 4] current state
                          beta,         # matrix [n_locs, 1] transmission
                          t_steps = 7,  # number of tau steps to simulate
                          tau     = 1,  # step size (days)
                          gamma   = 1/3,# I -> (H or R) rate
                          alpha   = 0.02, # fraction of I that go to H
                          eta     = 1/5) { # H -> R rate
  
  # Ensure `sihr` is always a matrix [n_locs, 4] even when n_locs = 1
  if (is.null(dim(sihr))) sihr <- matrix(sihr, nrow = 1)
  n_locs <- nrow(sihr)
  
  # Use a list of matrices rather than a 3D array to avoid R's dim-dropping
  # quirks when n_locs = 1. Each element is a [n_locs, 4] matrix.
  out <- vector("list", t_steps + 1)
  out[[1]] <- sihr
  
  state <- sihr
  for (step in seq_len(t_steps)) {
    S <- state[, 1]; I <- state[, 2]; H <- state[, 3]; R <- state[, 4]
    beta_vec <- as.numeric(beta)  # force to plain vector
    
    # Expected transition counts over tau
    lambda_SI <- pmax(beta_vec * S * I / population, 0) * tau
    lambda_IH <- pmax(gamma * alpha * I, 0) * tau
    lambda_IR <- pmax(gamma * (1 - alpha) * I, 0) * tau
    lambda_HR <- pmax(eta * H, 0) * tau
    
    dSI <- rpois(n_locs, lambda_SI)
    dIH <- rpois(n_locs, lambda_IH)
    dIR <- rpois(n_locs, lambda_IR)
    dHR <- rpois(n_locs, lambda_HR)
    
    # Cap so no compartment goes negative
    dSI <- pmin(dSI, S)
    total_I_out <- dIH + dIR
    over <- total_I_out > I
    if (any(over)) {
      avail <- I[over]
      tot   <- total_I_out[over]
      dIH[over] <- floor(dIH[over] * avail / tot)
      dIR[over] <- avail - dIH[over]
    }
    dHR <- pmin(dHR, H)
    
    new_state <- matrix(0, nrow = n_locs, ncol = 4)
    new_state[, 1] <- S - dSI
    new_state[, 2] <- I + dSI - dIH - dIR
    new_state[, 3] <- H + dIH - dHR
    new_state[, 4] <- R + dIR + dHR
    
    state <- new_state
    out[[step + 1]] <- state
  }
  
  # Convert list of [n_locs, 4] matrices to array [n_locs, 4, t_steps+1]
  result <- array(0, dim = c(n_locs, 4, t_steps + 1))
  for (s in seq_along(out)) {
    result[, , s] <- out[[s]]
  }
  result
}