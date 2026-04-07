## Modified from Galiana et al 2023 https://doi.org/10.1038/s41559-023-02254-y (data-analysis_occurrence.R)
# decided not to use this code

############ FITTING MODELS TO DEGREE DISTRIBUTIONS (Supplementary Table S2) ###############################

require(AICcmodavg)
require(dplyr)
require(purrr)
require(sars)
require(ggplot2)

####### example from stats::nls() documentation
# x <- -(1:100)/10
# y <- 100 + 10 * exp(x / 2) + rnorm(x)/10
# nlmod <- nls(y ~  Const + A * exp(B * x))

# plot(x,y, main = "nls(*), data, true function and fit, n=100")
# curve(100 + 10 * exp(x / 2), col = 4, add = TRUE)
# lines(x, predict(nlmod), col = 2)
######


# this is a cvs with columns called "species" and "interactions" where interactions is the node degree
output_galpar_real_pred <- read.csv("/home/fiona_callahan/TheoryCoOccur/network_real_pred_galpar.csv", header = TRUE)
output_galpar_real_pred <- output_galpar_real_pred[, c("species", "interactions")]

# distr = "rand" # "exp", "power", "rand"
# cutoff <- 0.1 # cutoffs = [0.05, 0.1, 0.15, 0.2]
# i=0

# fileNameCooc <- paste0("/home/fiona_callahan/TheoryCoOccur/degree_sequences/cooc/cooc_nodeDegree_", distr, "_cutoff", cutoff, "_", i, ".csv")
# fileNameInter <- paste0("/home/fiona_callahan/TheoryCoOccur/degree_sequences/inter/inter_nodeDegree_", distr, "_", i, ".csv")

# nodeDegree_df <- read.csv(fileNameCooc, header = TRUE)

df <- unique(output_galpar_real_pred) # Data corresponding to the network of biotic interactions
#df <- unique(nodeDegree_df) # Data corresponding to the network of biotic interactions
occur_pot_prey = as.vector(table(df$interactions)) # get counts of different numbers of interactions
occur_pot_prey = occur_pot_prey/sum(occur_pot_prey)
p = occur_pot_prey/sum(occur_pot_prey)
y = rev(cumsum(rev(p)))
x = as.numeric(names(table(df$interactions)))

temp_real <- data.frame(x, y)

## NETWORK OF BIOTIC INTERACTIONS

degree_dist_params <- NULL

failed <- tryCatch({
  mod1 <- tryCatch({ # power law with exponential cutoff
    nls(y ~ ( (x^-a) *(exp(-x/b))), data = temp_real, start = list(a = .0001, b = 2), control=nls.control(maxiter = 1e3))
  }, error = function(e) {
    NA
  }, finally = {
  })
  mod2 <- tryCatch({ # exponential
    nls(y ~ (exp(-x/b)), data = temp_real, start = list(b = 2), control=nls.control(maxiter = 1000))
  }, error = function(e) {
    NA
  }, finally = {
  })
  mod3 <- tryCatch({  # power law
    nls(y ~ (x^-a), data = temp_real, start = list(a = .01), control=nls.control(maxiter = 1e3))
  }, error = function(e) {
    NA
  }, finally = {
  })
  mod4 <- tryCatch({ # lognormal
    nls(y ~ ( (1/ (x * b * sqrt(2*pi) )) * exp(- ( ((log(x) - a)^2) / (2*(b^2)) )) ), data = temp_real, start = list(a = .3, b = .3), control=nls.control(maxiter = 1000))
  }, error = function(e) {
    NA
  }, finally = {
  })
  FALSE
}, warning = function(w) {
  FALSE
}, error = function(e) {
  TRUE
}, finally = {
  
})

if(!failed){
  model_list <- list(mod1,mod2,mod3,mod4)
  names_list <- c('mod1','mod2','mod3','mod4')
  if(length(which(is.na(model_list))) != 0){
    names_list <- names_list[-which(is.na(model_list))]
    model_list <- model_list[-which(is.na(model_list))]
  }
  names(model_list) <- names_list
  
  if(length(model_list) == 0) next
  
  if(length(model_list) == 1){
    model_name <- names_list[1]
    model <- summary(eval(as.symbol(model_name)))
  }else{
    aic_comp <- aictab(model_list)
    model_name <- tolower(as.character(aic_comp$Modnames[1]))
    model <- summary(eval(as.symbol(model_name)))
  }
  
  if(model_name == 'mod1' | model_name == 'mod4'){
    cur_out <- data.frame(model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=model$coefficients[2,1], b.std.err=model$coefficients[2,2], b.tval=model$coefficients[2,3], b.pval=model$coefficients[2,4])
  }else if(model_name == 'mod2'){
    cur_out <- data.frame( model=model_name, a=NA, a.std.err=NA, a.tval=NA, a.pval=NA, b=model$coefficients[1,1], b.std.err=model$coefficients[1,2], b.tval=model$coefficients[1,3], b.pval=model$coefficients[1,4])
  }else{
    cur_out <- data.frame(model=model_name, a=model$coefficients[1,1], a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3], a.pval=model$coefficients[1,4], b=NA, b.std.err=NA, b.tval=NA, b.pval=NA)
  }
  
  if(is.null(degree_dist_params)){
    degree_dist_params <- cur_out
  }else{
    degree_dist_params <- rbind(degree_dist_params, cur_out)
  }
}

aic_comp
degree_dist_params
