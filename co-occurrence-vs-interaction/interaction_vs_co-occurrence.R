library(ggplot2)

nspecies_gr1 <- 3
nspecies_gr2 <- 4

nested_mx <- matrix(c(
  1, 1, 0, 0,
  1, 1, 1, 0,
  1, 1, 1, 1
), nrow = nspecies_gr1, ncol = nspecies_gr2, byrow = TRUE)

# A here is from the model specidied in May 1972
A_from_bipartite <- function(bipartite_mx, type = "mutualist") {
    nspecies_gr1 <- dim(bipartite_mx)[1]
    nspecies_gr2 <- dim(bipartite_mx)[2]
    nSp_total <- nspecies_gr1 + nspecies_gr2
    A <- matrix(data = 0, nrow = nSp_total, ncol = nSp_total)
    # put interactions into matrix
    for (i in 1:nspecies_gr1) {
      for (j in 1:nspecies_gr2) {
        if (bipartite_mx[i, j] != 0) {
          if (type == "parasite") {
            A[i, nspecies_gr1 + j] <- runif(n = 1, 0, 1)
            A[nspecies_gr1 + j, i] <- runif(n = 1, -1, 0)
          } else if (type == "mutualist") {
            A[i, nspecies_gr1 + j] <- runif(n = 1, 0, 1)
            A[nspecies_gr1 + j, i] <- runif(n = 1, 0, 1)
          }
        }
      }
    }
    # diagonal = -1
    for (k in 1:nSp_total) {
      A[k, k] <- -1
    }
    A
}

run_steps <- function(A, x, n_steps = 1, verbose = FALSE) {
  for (i in 1:n_steps) {
    if (verbose) {
      print(x)
    }
    x <- x + A %*% x
    x <- x * (x > 0)
  }
   x
}

# simulating time steps plus noise to add disturbances to get a set of data points through "time"
run_simulation <- function(A, x_init, n = 100, burn = 100, sigma_noise = 0, skip_steps = 10) {
  n_species <- dim(A)[1]
  data <- data.frame(matrix(data = NA, nrow = n, ncol = 1 + n_species))
  names_sp <- sapply(1:n_species, FUN = function(x) {paste0("Sp", x)})
  names(data) <- c("t", names_sp)
  x <- run_steps(A = A, x = x_init, n_steps = burn)
  tt <- burn
  for (i in 1:n) {
    x <- x + rnorm(n = length(x), 0, sd = sigma_noise)
    x <- run_steps(A = A, x = x, n_steps = skip_steps)
    tt <- tt + skip_steps
    data[i, ] <- c(tt, x)
  }
  data
}


A <- A_from_bipartite(nested_mx, type = "mutualist")
A <- A_from_bipartite(nested_mx, type = "parasite")
eigen(A)$values
x <- matrix(data = runif(n = dim(A)[1], min = 0, max = 1000), nrow = dim(A)[1], ncol = 1)

run_steps(A, x, n_steps = 10, verbose = TRUE)

data1 <- run_simulation(A, x)



################# lotka voltera #######################
# all parms are positive
alpha <- 0.01 # effect of predator on prey
beta <- 0.005 # effect of prey on predator
a <- 1.5 # birth rate of prey
b <- 0.5 # death rate of predator

A <- matrix(data = c(0, -alpha * b / beta,
                    beta * a / alpha, 0), nrow = 2, ncol = 2, byrow = TRUE)

eigen(A)

# starting values
h_0 <- 51
p_0 <- 49

t_steps <- 100 # number of time steps
H <- rep(NA, times = t_steps) # prey abundance
P <- rep(NA, times = t_steps) # predator abundance

H[1] <- h_0
P[1] <- p_0

# I don't understand how this relates to the equations in May book p 41 
# I think it is wrong
# for (tt in 2:t_steps) {
#   change <- A %*% matrix(data = c(H[tt - 1], P[tt - 1]), ncol = 1)
#   H[tt] <- H[tt - 1] + change[1]
#   P[tt] <- P[tt - 1] + change[2]
# }
dt <- 0.0001
current_h <- h_0
current_p <- p_0
for (tt in 2:t_steps) {
  for (ttt in 1:round(1 / dt)) {
    dH <- current_h * (a - alpha * current_p)
    dP <- current_p * (-b + beta * current_h)
    current_h <- current_h + dH * dt
    current_p <- current_p + dP * dt
  }
  # store only occasionally
  H[tt] <- current_h
  P[tt] <- current_p
}

df <- data.frame(time = 1:t_steps, Prey = H, Predator = P)

# Plot time series
ggplot(df, aes(x = time)) +
  geom_line(aes(y = Prey, color = "Prey")) +
  geom_line(aes(y = Predator, color = "Predator")) +
  labs(title = "Lotka-Volterra Dynamics Over Time",
       y = "Population", color = "Species") +
  theme_minimal()

# Phase-plane plot (Predator vs Prey)
ggplot(df, aes(x = Prey, y = Predator)) +
  geom_path(color = "purple") +
  labs(title = "Phase Plane: Predator vs Prey",
       x = "Prey Population", y = "Predator Population") +
  theme_minimal()

