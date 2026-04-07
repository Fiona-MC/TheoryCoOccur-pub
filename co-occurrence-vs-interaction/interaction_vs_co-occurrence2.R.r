# Parameters
N <- 7           # Number of species
C <- 0.2          # Connectance
sigma <- 0.1      # Std dev of interaction strengths
d <- 1            # Self-regulation strength

# Generate random community matrix
A <- matrix(0, nrow = N, ncol = N)

# Fill off-diagonal elements with random interactions
for (i in 1:N) {
  for (j in 1:N) {
    if (i != j && runif(1) < C) {
      A[i, j] <- rnorm(1, mean = 0, sd = sigma)
    }
  }
}

# Set diagonal elements to -d (self-regulation)
diag(A) <- -d

# Eigenvalues
eig_vals <- eigen(A)$values

# Plot eigenvalues
plot(Re(eig_vals), Im(eig_vals), xlab = "Real part", ylab = "Imaginary part",
     main = "Eigenvalues of the Community Matrix", pch = 19, col = "blue")
abline(v = 0, col = "red", lty = 2)

A


