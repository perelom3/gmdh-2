library(purrr)
library(matrixcalc)
library(datastructures)
library(dplyr) # for the pipe


# ----INPUT FROM USER -----
# frequencies amount(m) and interval
interval_start <- 2
interval_end <- 11
# stop threshold 
d <- 0.0001
# D - ES epssilone
D <- 0.001
# GA operators
# start amount of applicants (M)
M <- 100
# k - interations to double check
k <- 5
# F - freedom paramu
F <- 5
# mutation probability
mutation_probability <- 0.2
# cross over probability
cross_over_probalitity <- 0.9
# seed
# H amount of model for next breed
H <- 80
# N - gen length & number of frequencies
N <- 50
Seed <- 123

# Beta coefficient
Beta <- 0.5
# Alpha coefficient
Alpha <- 0.5

# Read data
Data <<- read.csv(file="./data/test1.csv", header=FALSE, sep=",")


# Split to train test exam parts
train_size <- floor(0.70 * nrow(Data))
test_size <- floor(0.20 * nrow(Data))
exam_size <- floor(0.10 * nrow(Data))

set.seed(Seed)

train_ind <- sample(seq_len(nrow(Data)), size = train_size)
test_ind <- sample(seq_len(nrow(Data)), size = test_size)
exam_ind <- sample(seq_len(nrow(Data)), size = exam_size)

train <- Data[train_ind, ]
test <- Data[test_ind, ]
exam <-  Data[exam_ind, ]

# --- PREVIEW ----
# Data
# FFS over data

plot(Data, type="l")
Data.spec <- spectrum(Data[["V2"]],span=5,plot=FALSE)
spx <- Data.spec$freq * 44
spy <- 2*Data.spec$spec
plot(spy~spx,xlab="frequency",ylab="spectral density",type="l")

# ---- GMDH --------

X <- train[["V1"]]
Y <- train[["V2"]]


## ---- UTILS ----

print_list <- function(...) {
  print(c(...), collapse = " ")
}

create_gen <- function(N) {
  sample(c(0, 1), N ,replace=TRUE)
  base <- rep(0, N)
  base[sample(seq(1, N), 2, replace = TRUE)] <- 1
  base
}

lsm <- function (X, Y, formula) {
  #matrixcalc::matrix.inverse(t(X) %*% X) %*% t(X) %*% Y
  lm(formula)
}

make_X <- function( ws) {
    paste("Y", paste(unlist(purrr::map(ws, function(freq) {
      paste("sin(",freq," * X) + cos(",freq, "* X)", sep="")
    })), collapse = " + "), sep = "~")
}

make_model <- function(In, Out, gen) {
  ws <- get_freqs(fs, gen)
  formula <- make_X(ws)
  Y <- matrix(Out)
  X <- matrix(In)
  lsm_result <- lsm(X, Y, formula)
  coeff_l <- length(lsm_result$coefficients)
  model_coefficients <- unname(lsm_result$coefficients[2:coeff_l])
  constant_part <- lsm_result$coefficients[1]
  list(fs=ws, coeff=model_coefficients, gen=gen, lsm=lsm_result, const=constant_part)
}

get_freqs <- function(interval, gen) {
  interval[gen == 1]
}

compute_model_over <- function(data, model) {
  unlist(map(data, function(v) {
    sum(model$const, unlist(map(seq_along(model$fs), function(i) {
      model$coeff[(i - 1) * 2 + 1] * sin(model$fs[i] * v ) + model$coeff[(i - 1) * 2 + 2] * cos(model$fs[i] * v )
    })))
  }))
}

compute_delta_over <- function(real_out, model_out, real_mean) {
  distance <- sum(unlist(map2(real_out, model_out, function(r, o) {
    (o - r) ^ 2 
  })))
  
  std <- sum(unlist(map(model_out, function(o) {
    (o - real_mean) ^ 2 
  })))
  
  sqrt(distance / std)
}

compute_ec <- function(model, train_in, train_out, test_in, test_out, Aplha, Beta) {
  delta_A_A <- compute_delta_over(train_out, compute_model_over(train_in, model), mean(train_out))
  delta_A_B <- compute_delta_over(test_out, compute_model_over(test_in, model), mean(test_out))
  K <- Beta * (Alpha * delta_A_A + (1 - Alpha) * delta_A_B) + (1 - Beta) * (abs(delta_A_A - delta_A_B) / abs(delta_A_A + delta_A_B))
  model$K <- K
  model
}

simple_crossover <- function(left, right) {
  l <- max(length(left), length(right))
  split_point <- floor(runif(1, 2, l))
  
  first <- c(left[1:split_point], right[(split_point + 1):l])
  second <- c(right[1:split_point], left[(split_point + 1):l])
  list(first, second)
}

mutate <- function(gen) {
  indexs <- sample(seq_along(gen), 3, replace = TRUE)
  to_mutate <- gen[indexs]
  gen[indexs] <- !to_mutate
  list(gen)
}

breed_models <- function(selected_models_for_breed, mutation_probability, cross_over_probalitity, M) {
  breed <- list()
  
  for(i in 1:M) {
    choice <- runif(1, 0, 1)
    if (choice > (1 - cross_over_probalitity)) {
      breed <- append(breed, simple_crossover(sample(selected_models_for_breed, 1)[[1]]$gen, sample(selected_models_for_breed, 1)[[1]]$gen))
    }
    
    if (choice > (1 - mutation_probability)) {
      breed <- append(breed, mutate(sample(selected_models_for_breed, 1)[[1]]$gen))
    }
  }
  
  breed
}

#

saved_models <- datastructures::binomial_heap("numeric")

k_c <- k
prev_ca <- 0.0

# Generate frequenscies interval
fs <- unlist(map(seq(1, N), function(i) {
  ((interval_end  * i) / N)  + runif(1, interval_start /  N, interval_end / N)
})) 

iteration <- 0
# Generate M applicants with m chromosomes

gens = map(seq(0, M),function(x) { create_gen(N) })

ga_tick <- function(gens) {
  # -- compute least squares method(LSM)
  models <- map(gens, function(gen) {make_model(X, Y, gen)})
  # -- compute external criterion(EC) for model on test data
  models_with_ec <- map(models, function(model) {
    compute_ec(model, train[["V1"]], train[["V2"]], test[["V1"]], test[["V2"]], Aplha, Beta ) 
  })
  
  models_sorted_by_es <- models_with_ec[order(sapply(models_with_ec, function(x) x$K))]
  # select H < M best applicants for ne xt breed by EC
  selected_models_for_breed <- models_sorted_by_es[1:H]
  # Save F models
  map(selected_models_for_breed[1:F], function(model) {
    datastructures::insert(saved_models, model$K, model)
  })
  print_list(" [models] ", length(models))
  # breed new M > H applacants from H selected
  gens <- breed_models(selected_models_for_breed, mutation_probability, cross_over_probalitity, M)
  # Check if GA should stop by absolute value of difference between current best EC and previous, but keep going for k steps to double check
  list(K=selected_models_for_breed[[1]]$K, gens=gens)
}

while (TRUE) {
  # for each applicant
  # remove all zeros gens
  gens <<- gens[lapply(gens, sum) > 0]
  ga <- ga_tick(gens)
  current_ca <- ga$K
  gens <<- ga$gens
  
  print_list("Current minimum CA: ", current_ca, " [Iteration] ", iteration)

  if (abs(prev_ca - current_ca) < d) {
    k_c <<- k_c - 1
  } else {
    k_c <<- k
  }

  if (k_c == 0) {
    break
  }
  
  prev_ca <<- current_ca
  iteration <<- iteration + 1
}

# ---- Result ------
# print best 5 models
# plot best model and data
best_model <- datastructures::peek(saved_models)[[1]]
DD <- data.frame(model=compute_model_over(Data[["V1"]], best_model), real=Data[["V2"]])
matplot(DD, type="l")
plot(Data)
plot( data.frame(Data[["V1"]], compute_model_over(Data[["V1"]], best_model)))
# calc criterions on exam data
  