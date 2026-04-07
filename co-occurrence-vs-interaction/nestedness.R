# nestedness in interaction matrix =? nestedness in correlation matrix?
library(bipartite) # for nestedness metrics
library(igraph)
library(data.table) # for fast reading

# Here are the simulation matrices
nested_mx <- matrix(c(
    1, 0, 0, 0, 0,
    1, 1, 0, 0, 0,
    1, 1, 1, 0, 0,
    1, 1, 1, 1, 0,
    1, 1, 1, 1, 1
), nrow = 5, ncol = 5, byrow = TRUE)

modular_mx <- matrix(c(
    1, 1, 1, 0, 0,
    1, 1, 1, 0, 0,
    1, 1, 1, 0, 0,
    0, 0, 0, 1, 1,
    0, 0, 0, 1, 1
), nrow = 5, ncol = 5, byrow = TRUE)

alpha_from_bipartite <- function(bipartite_mx, type = "mutualist") {
    nspecies_gr1 <- dim(bipartite_mx)[1]
    nspecies_gr2 <- dim(bipartite_mx)[2]
    nSp_total <- nspecies_gr1 + nspecies_gr2
    alpha <- matrix(data = 0, nrow = nSp_total, ncol = nSp_total)
    # put interactions into matrix
    for (i in 1:nspecies_gr1) {
      for (j in 1:nspecies_gr2) {
        if (bipartite_mx[i, j] != 0) {
          if (type == "parasite") {
            alpha[i, nspecies_gr1 + j] <- 1
            alpha[nspecies_gr1 + j, i] <- -1
          } else if (type == "mutualist") {
            alpha[i, nspecies_gr1 + j] <- 1
            alpha[nspecies_gr1 + j, i] <- 1
          }
        }
      }
    }
    alpha
}

bipartite_from_alpha <- function(alpha, split_index) {
    nspecies_gr1 <- split_index
    nspecies_gr2 <- dim(alpha)[1] - split_index
    bipartite <- matrix(data = 0, nrow = nspecies_gr1, ncol = nspecies_gr2)
    for (i in 1:nspecies_gr1) {
      for (j in 1:nspecies_gr2) {
        bipartite[i, j] <- alpha[i, nspecies_gr1 + j]
      }
    }
    bipartite
}

nSpecies <- 10
sp_names <- sapply(1:nSpecies, function(x) {paste0("Sp", x)})

threshold <- 0.1
nested_res <- data.frame()
for (run in 1:5) {
  data <- as.data.frame(fread(file = paste0("/space/s1/fiona_callahan/networks_testSim/nested/mutualist/randomRun", run, "/sim_sitetab_readAbd_sampled.csv")))
  data_sp <- data[, sp_names]
  correlation_matrix <- cor(data_sp, method = "pearson")

  bipartite_mx <- bipartite_from_alpha(correlation_matrix, 5)
  bipartite_mx_binary <- (bipartite_mx > threshold) * 1
  bipartite_mx_binary <- bipartite_mx_binary[apply(bipartite_mx_binary, FUN = sum, MARGIN = 1) > 0,
                                              apply(bipartite_mx_binary, FUN = sum, MARGIN = 2) > 0]

  # NODF -- 0 to 100, 100 is very nested
  nested_res <- rbind(nested_res, nested(bipartite_mx_binary, method = "ALL", rescale = TRUE))
  names(nested_res) <- names(nested(bipartite_mx_binary, method = "ALL", rescale = TRUE))
}

mod_res <- data.frame()
for (run in 1:5) {
  data_modular <- as.data.frame(fread(file = paste0("/space/s1/fiona_callahan/networks_testSim/modular/mutualist/randomRun", run, "/sim_sitetab_readAbd_sampled.csv")))
  data_sp_modular <- data_modular[, sp_names]
  correlation_matrix_modular <- cor(data_sp_modular, method = "pearson")

  bipartite_mx_modular <- bipartite_from_alpha(correlation_matrix_modular, 5)
  bipartite_mx_mod_binary <- bipartite_mx_modular > threshold
  bipartite_mx_mod_binary <- bipartite_mx_mod_binary[apply(bipartite_mx_mod_binary, FUN = sum, MARGIN = 1) > 0,
                                              apply(bipartite_mx_mod_binary, FUN = sum, MARGIN = 2) > 0]

  mod_res <- rbind(mod_res, nested(bipartite_mx_mod_binary, method = "ALL", rescale = TRUE))
  names(mod_res) <- names(nested(bipartite_mx_mod_binary, method = "ALL", rescale = TRUE))
}

nested_res
mod_res

######################################################
### another little scratch analysis not sim-data-based
######################################################
adjacency_mx <- alpha_from_bipartite(nested_mx, type = "mutualist")
g <- graph_from_adjacency_matrix(adjacency_mx, mode = "undirected")
shortest_paths_matrix <- distances(g)

corr <- 0.5^shortest_paths_matrix
bipartite_mx_corr <- bipartite_from_alpha(corr, 5)
nested(bipartite_mx_corr, method = "NODF2")


adjacency_mx_mod <- alpha_from_bipartite(modular_mx, type = "mutualist")
g_mod <- graph_from_adjacency_matrix(adjacency_mx_mod, mode = "undirected")
shortest_paths_matrix_mod <- distances(g_mod)

corr_mod <- 0.5^shortest_paths_matrix_mod
bipartite_mx_mod <- bipartite_from_alpha(corr_mod, 5)
nested(bipartite_mx_mod, method = "NODF2")
