library(doRNG)
set.seed(1)
# Number of modes
d <- 3
# Number of users
n1 <- 100
# Number of contexts
n2 <- 9
# Number of items
n3 <- 100
# Number of time points
tt <- (tt1 <- 12) + (tt2 <- 8)
# Number of subgroups for users
m1 <- 10
# Number of subgroups for contexts
m2 <- 3
# Number of subgroups for items
m3 <- 10
# Number of subgroups for time
m4 <- 4
# Number of latent factors
r <- 3

s1 <- sample.int(m1, n1, replace = TRUE)
s2 <- sample.int(m2, n2, replace = TRUE)
s3 <- sample.int(m3, n3, replace = TRUE)
s4 <- sample.int(m4, tt, replace = TRUE)

t_horizon <- sort(runif(tt, 0, 1))

p1 <- mvtnorm::rmvnorm(n1, rep(0, r), diag(1, r))
p2 <- mvtnorm::rmvnorm(n2, rep(0, r), diag(1, r))
p3 <- mvtnorm::rmvnorm(n3, rep(0, r), diag(1, r))

h <- rep(list(NA), r)
h[[1]] <- function(t) {
  sin(0.3 * pi * t)
}
h[[2]] <- function(t) {
  8 * t * (1 - t) - 1
}
h[[3]] <- function(t) {
  cos(0.2 * pi * t) + 1
}
g <- rep(list(NA), m4)
g[[1]] <- function(t) {
  2 * t - 1
}
g[[2]] <- function(t) {
  8 * (t - 0.5) ^ 3
}
g[[3]] <- function(t) {
  sin(0.1 * pi * t) + cos(pi * t)
}
g[[4]] <- function(t) {
  -5 * exp(t) + 10
}
gfun <- function(t) {
  g[[s4[t == t_horizon]]](t)
}

q1 <- -1 + 0.4 * seq_len(m1)
q2 <- -1.2 + 0.6 * seq_len(m2)
q3 <- -0.4 + 0.2 * seq_len(m3)

epsilon <- rnorm(n1 * n2 * n3 * tt, 0, 1)
epsilon <- array(epsilon, c(n1, n2, n3, tt))

# y <- array(NA, c(n1, n2, n3, tt))
# for (i1 in seq_len(nrow(p1))) {
#   for (i2 in seq_len(nrow(p2))) {
#     for (i3 in seq_len(nrow(p3))) {
#       for (t in seq_len(tt)) {
#         y[i1, i2, i3, t] <- sum(
#           sapply(seq_len(r), function(j) {
#             h[[j]](t_horizon[t]) * p1[i1, j] * p2[i2, j] * p3[i3, j]
#           })
#         ) + gfun(t_horizon[t]) * q1[s1[i1]] * q2[s2[i2]] * q3[s3[i3]] + epsilon[i1, i2, i3, t]
#       }
#     }
#   }
# }

hmat <- matrix(NA, tt, r)
for (t in seq_len(tt)) {
  for (j in seq_len(r)) {
    hmat[t, j] <- h[[j]](t_horizon[t])
  }
}

decomp <- array(0, c(n1, n2, n3, r))
for (j in seq_len(r)) {
  decomp[, , , j] <- p1[, j] %o% p2[, j] %o% p3[, j]
}
qtensor <- q1[s1] %o% q2[s2] %o% q3[s3]

y <- array(0, c(n1, n2, n3, tt))
for (t in seq_len(tt)) {
  for (j in seq_len(r)) {
    y[, , , t] <- y[, , , t] + hmat[t, j] * decomp[, , , j]
  }
  y[, , , t] <- y[, , , t] + gfun(t_horizon[t]) * qtensor + epsilon[, , , t]
}

missing_rate <- 0.8
missing_tensor <- array(rbinom(n1 * n2 * n3 * tt, 1, 1 - missing_rate), c(n1, n2, n3, tt))

y <- y * missing_tensor

# Number of knots
an <- 3
# Continuous derivatives
kappa <- 2
nu <- NULL
omega <- NULL

mm <- an + kappa + 1

bb <- rep(list(rep(list(NA), mm)), r)
aa <- rep(list(rep(list(NA), mm)), m4)

for (j in seq_len(r)) {
  for (i in seq_len(mm)) {
    bb[[j]][[i]] <- function(t, j, i, kappa, nu) {
      if (i <= kappa + 1) {
        t^(i - 1)
      } else {
        max(t - nu[j, i - kappa - 1], 0)^kappa
      }
    }
  }
}

for (j in seq_len(m4)) {
  for (i in seq_len(mm)) {
    aa[[j]][[i]] <- function(t, j, i, kappa, omega) {
      if (i <= kappa + 1) {
        t^(i - 1)
      } else {
        max(t - omega[j, i - kappa - 1], 0)^kappa
      }
    }
  }
}


alpha <- matrix(0, r, mm)
beta <- matrix(0, m4, mm)

hhat <- rep(list(NA), r)
for (j in seq_len(r)) {
  hhat[[j]] <- function(t, j, kappa, alpha, bb, nu, mm) {
    # print(mm)
    sum(sapply(seq_len(mm), function(i) {
      alpha[j, i] * bb[[j]][[i]](t, j, i, kappa, nu)
    }))
  }
}

ghat <- rep(list(NA), m4)
for (j in seq_len(m4)) {
  ghat[[j]] <- function(t, j, kappa, beta, aa, omega, mm) {
    sum(sapply(seq_len(mm), function(i) {
      beta[j, i] * aa[[j]][[i]](t, j, i, kappa, omega)
    }))
  }
}
# ghatfun <- function(t, kappa, beta) {
#   ghat[[s4[t == t_horizon]]](t, kappa, beta)
# }

likelihood <- function(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) {
  yhat <- array(0, c(nrow(p1), nrow(p2), nrow(p3), tt))
  for (t in seq_len(tt)) {
    for (j in seq_len(r)) {
      yhat[, , , t] <- yhat[, , , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j] %o% p3[, j]
    }
    yhat[, , , t] <- yhat[, , , t] + ghat[[s4[t]]](t, s4[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2] %o% q3[s3]
  }
  indicator <- y != 0
  omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[4]), function(t) indicator[, , , t]))
  ll <- 0
  for (i1 in seq_len(nrow(p1))) {
    for (i2 in seq_len(nrow(p2))) {
      for (i3 in seq_len(nrow(p3))) {
        if (omega_set[i1, i2, i3]) {
          ll <- ll + (y[i1, i2, i3, ] - yhat[i1, i2, i3, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, i3, ] - yhat[i1, i2, i3, ])
        }
      }
    }
  }
  ll <- ll + lambda * (norm(p1, "F")^2 + norm(p2, "F")^2 + norm(p3, "F")^2 + norm(q1[s1], "2")^2 + norm(q2[s2], "2")^2 + norm(q3[s3], "2")^2 + norm(alpha, "F")^2 + norm(beta, "F")^2)
  ll
}

dtrs <- function(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, epsilon = 1e-4) {
  iteration_l <- 0
  error <- 1
  while (error >= epsilon) {
    iteration_l <- iteration_l + 1
    print(paste("Iteration", iteration_l))
    # Estimate p1
    print("Estimate p1")
    cl <- parallel::makeCluster(parallel::detectCores(), outfile = "log/p1.log")
    doParallel::registerDoParallel(cl)
    i1 <- NULL
    p1new <- doRNG::`%dorng%`(foreach::foreach(i1 = seq_len(nrow(p1)), .combine = "rbind"), {
      optim(
        par = p1[i1, ],
        fn = function(p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, r, y, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
          yhat <- array(0, c(nrow(p2), nrow(p3), length(t_horizon)))
          for (t in seq_len(length(t_horizon))) {
            for (j in seq_len(r)) {
              yhat[, , t] <- yhat[, , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[j] * p2[, j] %o% p3[, j]
            }
            yhat[, , t] <- yhat[, , t] + ghat[[s4[t]]](t, s4[t], kappa, beta, aa, omega, mm) * q1[s1][i1] * q2[s2] %o% q3[s3]
          }
          indicator <- y != 0
          omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[3]), function(x) indicator[ , , x]))
          result <- 0
          for (i2 in seq_len(nrow(p2))) {
            for (i3 in seq_len(nrow(p3))) {
              if (omega_set[i2, i3]) {
                result <- result + (y[i2, i3, ] - yhat[i2, i3, ]) %*% diag(length(t_horizon)) %*% (y[i2, i3, ] - yhat[i2, i3, ])
              }
            }
          }
          result <- result + lambda * norm(p1, type = "2")^2
        },
        p2 = p2,
        p3 = p3,
        q1 = q1,
        q2 = q2,
        q3 = q3,
        s1 = s1,
        s2 = s2,
        s3 = s3,
        s4 = s4,
        t_horizon = t_horizon,
        r = r,
        y = y[i1, , , ],
        alpha = alpha,
        beta = beta,
        hhat = hhat,
        ghat = ghat,
        aa = aa,
        bb = bb,
        mm = mm,
        nu = nu,
        omega = omega,
        lambda = lambda
      )$par
    })
    parallel::stopCluster(cl)
    jp1 <- 1 - likelihood(y, r, kappa, lambda, p1new, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)
    # Estimate p2
    print("Estimate p2")
    cl <- parallel::makeCluster(parallel::detectCores(), outfile = "log/p2.log")
    doParallel::registerDoParallel(cl)
    i2 <- NULL
    p2new <- doRNG::`%dorng%`(foreach::foreach(i2 = seq_len(nrow(p2)), .combine = "rbind"), {
      optim(
        par = p2[i2, ],
        fn = function(p2, p1, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, r, y, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
          yhat <- array(0, c(nrow(p1), nrow(p3), length(t_horizon)))
          for (t in seq_len(length(t_horizon))) {
            for (j in seq_len(r)) {
              yhat[, , t] <- yhat[, , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p2[j] * p1[, j] %o% p3[, j]
            }
            yhat[, , t] <- yhat[, , t] + ghat[[s4[t]]](t, s4[t], kappa, beta, aa, omega, mm) * q2[s2][i2] * q1[s1] %o% q3[s3]
          }
          indicator <- y != 0
          omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[3]), function(x) indicator[ , , x]))
          result <- 0
          for (i1 in seq_len(nrow(p1))) {
            for (i3 in seq_len(nrow(p3))) {
              if (omega_set[i1, i3]) {
                result <- result + (y[i1, i3, ] - yhat[i1, i3, ]) %*% diag(length(t_horizon)) %*% (y[i1, i3, ] - yhat[i1, i3, ])
              }
            }
          }
          result <- result + lambda * norm(p2, type = "2")^2
        },
        p1 = p1,
        p3 = p3,
        q1 = q1,
        q2 = q2,
        q3 = q3,
        s1 = s1,
        s2 = s2,
        s3 = s3,
        s4 = s4,
        t_horizon = t_horizon,
        r = r,
        y = y[, i2, , ],
        alpha = alpha,
        beta = beta,
        hhat = hhat,
        ghat = ghat,
        aa = aa,
        bb = bb,
        mm = mm,
        nu = nu,
        omega = omega,
        lambda = lambda
      )$par
    })
    parallel::stopCluster(cl)
    jp2 <- 1 - likelihood(y, r, kappa, lambda, p1, p2new, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)
    # Estimate p3
    print("Estimate p3")
    cl <- parallel::makeCluster(parallel::detectCores(), outfile = "log/p3.log")
    doParallel::registerDoParallel(cl)
    i3 <- NULL
    p3new <- doRNG::`%dorng%`(foreach::foreach(i3 = seq_len(nrow(p3)), .combine = "rbind"), {
      optim(
        par = p3[i3, ],
        fn = function(p3, p1, p2, q1, q2, q3, s1, s2, s3, s4, t_horizon, r, y, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
          yhat <- array(0, c(nrow(p1), nrow(p2), length(t_horizon)))
          for (t in seq_len(length(t_horizon))) {
            for (j in seq_len(r)) {
              yhat[, , t] <- yhat[, , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p3[j] * p1[, j] %o% p2[, j]
            }
            yhat[, , t] <- yhat[, , t] + ghat[[s4[t]]](t, s4[t], kappa, beta, aa, omega, mm) * q3[s3][i3] * q1[s1] %o% q2[s2]
          }
          indicator <- y != 0
          omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[3]), function(x) indicator[ , , x]))
          result <- 0
          for (i1 in seq_len(nrow(p1))) {
            for (i2 in seq_len(nrow(p2))) {
              if (omega_set[i1, i2]) {
                result <- result + (y[i1, i2, ] - yhat[i1, i2, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, ] - yhat[i1, i2, ])
              }
            }
          }
          result <- result + lambda * norm(p2, type = "2")^2
        },
        p1 = p1,
        p2 = p2,
        q1 = q1,
        q2 = q2,
        q3 = q3,
        s1 = s1,
        s2 = s2,
        s3 = s3,
        s4 = s4,
        t_horizon = t_horizon,
        r = r,
        y = y[, , i3, ],
        alpha = alpha,
        beta = beta,
        hhat = hhat,
        ghat = ghat,
        aa = aa,
        bb = bb,
        mm = mm,
        nu = nu,
        omega = omega,
        lambda = lambda
      )$par
    })
    parallel::stopCluster(cl)
    jp3 <- 1 - likelihood(y, r, kappa, lambda, p1, p2, p3new, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)
    # Estimate alpha
    print("Estimate alpha")
    alphanew <- optim(
      par = c(alpha),
      fn = function(alpha, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, r, y, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
        alpha <- matrix(alpha, ncol = mm)
        yhat <- array(0, c(nrow(p1), nrow(p2), nrow(p3), length(t_horizon)))
        for (t in seq_len(length(t_horizon))) {
          for (j in seq_len(r)) {
            yhat[, , , t] <- yhat[, , , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j] %o% p3[, j]
          }
          yhat[, , , t] <- yhat[, , , t] + ghat[[s4[t]]](t, s4[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2] %o% q3[s3]
        }
        indicator <- y != 0
        omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[4]), function(x) indicator[ , , , x]))
        result <- 0
        for (i1 in seq_len(nrow(p1))) {
          for (i2 in seq_len(nrow(p2))) {
            for (i3 in seq_len(nrow(p3))) {
              if (omega_set[i1, i2, i3]) {
                result <- result + (y[i1, i2, i3, ] - yhat[i1, i2, i3, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, i3, ] - yhat[i1, i2, i3, ])
              }
            }
          }
        }
        result <- result + lambda * norm(alpha, type = "F")^2
      },
      p1 = p1,
      p2 = p2,
      p3 = p3,
      q1 = q1,
      q2 = q2,
      q3 = q3,
      s1 = s1,
      s2 = s2,
      s3 = s3,
      s4 = s4,
      t_horizon = t_horizon,
      r = r,
      y = y,
      beta = beta,
      hhat = hhat,
      ghat = ghat,
      aa = aa,
      bb = bb,
      mm = mm,
      nu = nu,
      omega = omega,
      lambda = lambda
    )$par
    alphanew <- matrix(alphanew, ncol = mm)
    jalpha <- 1 - likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alphanew, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)

    jmax_index <- which.max(c(jp1, jp2, jp3, jalpha))
    if (jmax_index == 1) {
      p1 <- p1new
    } else if (jmax_index == 2) {
      p2 <- p2new
    } else if (jmax_index == 3) {
      p3 <- p3new
    } else if (jmax_index == 4) {
      alpha <- alphanew
    }
    # Estimate q1
    print("Estimate q1")
    q1new <- optim(
      par = q1,
      fn = function(q1, p1, p2, p3, q2, q3, s1, s2, s3, s4, t_horizon, r, y, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
        yhat <- array(0, c(nrow(p1), nrow(p2), nrow(p3), length(t_horizon)))
        for (t in seq_len(length(t_horizon))) {
          for (j in seq_len(r)) {
            yhat[, , , t] <- yhat[, , , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j] %o% p3[, j]
          }
          yhat[, , , t] <- yhat[, , , t] + ghat[[s4[t]]](t, s4[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2] %o% q3[s3]
        }
        indicator <- y != 0
        omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[4]), function(x) indicator[ , , , x]))
        result <- 0
        for (i1 in seq_len(nrow(p1))) {
          for (i2 in seq_len(nrow(p2))) {
            for (i3 in seq_len(nrow(p3))) {
              if (omega_set[i1, i2, i3]) {
                result <- result + (y[i1, i2, i3, ] - yhat[i1, i2, i3, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, i3, ] - yhat[i1, i2, i3, ])
              }
            }
          }
        }
        result <- result + lambda * norm(q1, type = "2")^2
      },
      p1 = p1,
      p2 = p2,
      p3 = p3,
      q2 = q2,
      q3 = q3,
      s1 = s1,
      s2 = s2,
      s3 = s3,
      s4 = s4,
      t_horizon = t_horizon,
      r = r,
      y = y,
      alpha = alpha,
      beta = beta,
      hhat = hhat,
      ghat = ghat,
      aa = aa,
      bb = bb,
      mm = mm,
      nu = nu,
      omega = omega,
      lambda = lambda
    )$par
    jq1 <- 1 - likelihood(y, r, kappa, lambda, p1, p2, p3, q1new, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)
    # Estimate q2
    print("Estimate q2")
    q2new <- optim(
      par = q2,
      fn = function(q2, p1, p2, p3, q1, q3, s1, s2, s3, s4, t_horizon, r, y, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
        yhat <- array(0, c(nrow(p1), nrow(p2), nrow(p3), length(t_horizon)))
        for (t in seq_len(length(t_horizon))) {
          for (j in seq_len(r)) {
            yhat[, , , t] <- yhat[, , , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j] %o% p3[, j]
          }
          yhat[, , , t] <- yhat[, , , t] + ghat[[s4[t]]](t, s4[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2] %o% q3[s3]
        }
        indicator <- y != 0
        omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[4]), function(x) indicator[ , , , x]))
        result <- 0
        for (i1 in seq_len(nrow(p1))) {
          for (i2 in seq_len(nrow(p2))) {
            for (i3 in seq_len(nrow(p3))) {
              if (omega_set[i1, i2, i3]) {
                result <- result + (y[i1, i2, i3, ] - yhat[i1, i2, i3, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, i3, ] - yhat[i1, i2, i3, ])
              }
            }
          }
        }
        result <- result + lambda * norm(q2, type = "2")^2
      },
      p1 = p1,
      p2 = p2,
      p3 = p3,
      q1 = q1,
      q3 = q3,
      s1 = s1,
      s2 = s2,
      s3 = s3,
      s4 = s4,
      t_horizon = t_horizon,
      r = r,
      y = y,
      alpha = alpha,
      beta = beta,
      hhat = hhat,
      ghat = ghat,
      aa = aa,
      bb = bb,
      mm = mm,
      nu = nu,
      omega = omega,
      lambda = lambda
    )$par
    jq2 <- 1 - likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2new, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)
    # Estimate q3
    print("Estimate q3")
    q3new <- optim(
      par = q3,
      fn = function(q3, p1, p2, p3, q1, q2, s1, s2, s3, s4, t_horizon, r, y, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
        yhat <- array(0, c(nrow(p1), nrow(p2), nrow(p3), length(t_horizon)))
        for (t in seq_len(length(t_horizon))) {
          for (j in seq_len(r)) {
            yhat[, , , t] <- yhat[, , , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j] %o% p3[, j]
          }
          yhat[, , , t] <- yhat[, , , t] + ghat[[s4[t]]](t, s4[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2] %o% q3[s3]
        }
        indicator <- y != 0
        omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[4]), function(x) indicator[ , , , x]))
        result <- 0
        for (i1 in seq_len(nrow(p1))) {
          for (i2 in seq_len(nrow(p2))) {
            for (i3 in seq_len(nrow(p3))) {
              if (omega_set[i1, i2, i3]) {
                result <- result + (y[i1, i2, i3, ] - yhat[i1, i2, i3, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, i3, ] - yhat[i1, i2, i3, ])
              }
            }
          }
        }
        result <- result + lambda * norm(q3, type = "2")^2
      },
      p1 = p1,
      p2 = p2,
      p3 = p3,
      q1 = q1,
      q2 = q2,
      s1 = s1,
      s2 = s2,
      s3 = s3,
      s4 = s4,
      t_horizon = t_horizon,
      r = r,
      y = y,
      alpha = alpha,
      beta = beta,
      hhat = hhat,
      ghat = ghat,
      aa = aa,
      bb = bb,
      mm = mm,
      nu = nu,
      omega = omega,
      lambda = lambda
    )$par
    jq3 <- 1 - likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3new, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)
    # Estimate beta
    print("Estimate beta")
    betanew <- optim(
      par = c(beta),
      fn = function(beta, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, r, y, alpha, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
        beta <- matrix(beta, ncol = mm)
        yhat <- array(0, c(nrow(p1), nrow(p2), nrow(p3), length(t_horizon)))
        for (t in seq_len(length(t_horizon))) {
          for (j in seq_len(r)) {
            yhat[, , , t] <- yhat[, , , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j] %o% p3[, j]
          }
          yhat[, , , t] <- yhat[, , , t] + ghat[[s4[t]]](t, s4[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2] %o% q3[s3]
        }
        indicator <- y != 0
        omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[4]), function(x) indicator[ , , , x]))
        result <- 0
        for (i1 in seq_len(nrow(p1))) {
          for (i2 in seq_len(nrow(p2))) {
            for (i3 in seq_len(nrow(p3))) {
              if (omega_set[i1, i2, i3]) {
                result <- result + (y[i1, i2, i3, ] - yhat[i1, i2, i3, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, i3, ] - yhat[i1, i2, i3, ])
              }
            }
          }
        }
        result <- result + lambda * norm(beta, type = "2")^2
      },
      p1 = p1,
      p2 = p2,
      p3 = p3,
      q1 = q1,
      q2 = q2,
      q3 = q3,
      s1 = s1,
      s2 = s2,
      s3 = s3,
      s4 = s4,
      t_horizon = t_horizon,
      r = r,
      y = y,
      alpha = alpha,
      hhat = hhat,
      ghat = ghat,
      aa = aa,
      bb = bb,
      mm = mm,
      nu = nu,
      omega = omega,
      lambda = lambda
    )$par
    betanew <- matrix(betanew, ncol = mm)
    jbeta <- 1 - likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, betanew, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)

    jmax_index <- which.max(c(jq1, jq2, jq3, jbeta))
    if (jmax_index == 1) {
      q1 <- q1new
    } else if (jmax_index == 2) {
      q2 <- q2new
    } else if (jmax_index == 3) {
      q3 <- q3new
    } else if (jmax_index == 4) {
      beta <- betanew
    }
    error <- max(jp1, jp2, jp3, jalpha, jq1, jq2, jq3, jbeta)
    print(c(jp1, jp2, jp3, jalpha, jq1, jq2, jq3, jbeta))
    saveRDS(list(p1 = p1, p2 = p2, p3 = p3, q1 = q1, q2 = q2, q3 = q3, alpha = alpha, beta = beta), paste0("dtrs", iteration_l, ".rds"))
  }
  return(list(p1 = p1, p2 = p2, p3 = p3, q1 = q1, q2 = q2, q3 = q3, alpha = alpha, beta = beta))
}

# Assmue to be known: r, s1, s2, s3, s4
p1 <- matrix(rnorm(n1 * r), n1, r)
p2 <- matrix(rnorm(n2 * r), n2, r)
p3 <- matrix(rnorm(n3 * r), n3, r)
q1 <- rnorm(m1)
q2 <- rnorm(m2)
q3 <- rnorm(m3)
alpha <- matrix(0, r, mm)
beta <- matrix(0, m4, mm)

result <- dtrs(y, r, kappa, lambda = 5, p1, p2, p3, q1, q2, q3, s1, s2, s3, s4, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, epsilon = 1e-4)

#> [1] "Iteration 1"
#> [1] "Estimate p1"
#> [1] "Estimate p2"
#> [1] "Estimate p3"
#> [1] "Estimate alpha"
#> [1] "Estimate q1"
#> [1] "Estimate q2"
#> [1] "Estimate q3"
#> [1] "Estimate beta"
#> [1] 5.968708e-04 6.862882e-05 0.000000e+00 2.475967e-06 2.847993e-04 3.324999e-05 1.644449e-04 1.614613e-04
#> [1] "Iteration 2"
#> [1] "Estimate p1"
#> [1] "Estimate p2"
#> [1] "Estimate p3"
#> [1] "Estimate alpha"
#> [1] "Estimate q1"
#> [1] "Estimate q2"
#> [1] "Estimate q3"
#> [1] "Estimate beta"
#> [1] 9.637809e-10 6.868937e-05 0.000000e+00 5.284007e-06 4.009411e-06 3.326175e-05 1.645031e-04 2.613922e-04
#> [1] "Iteration 3"
#> [1] "Estimate p1"
#> [1] "Estimate p2"
#> [1] "Estimate p3"
#> [1] "Estimate alpha"
#> [1] "Estimate q1"
#> [1] "Estimate q2"
#> [1] "Estimate q3"
#> [1] "Estimate beta"
#> [1] 9.640991e-10 2.962269e-09 0.000000e+00 2.251005e-09 1.489483e-04 8.562161e-04 1.164403e-03 2.407839e-06
#> [1] "Iteration 4"
#> [1] "Estimate p1"
#> [1] "Estimate p2"
#> [1] "Estimate p3"
#> [1] "Estimate alpha"
#> [1] "Estimate q1"
#> [1] "Estimate q2"
#> [1] "Estimate q3"
#> [1] "Estimate beta"
#> [1] 9.652230e-10 2.405397e-09 0.000000e+00 2.137092e-09 1.552279e-03 1.192340e-02 2.223632e-05 6.950460e-05
#> [1] "Iteration 5"
#> [1] "Estimate p1"
#> [1] "Estimate p2"
#> [1] "Estimate p3"
#> [1] "Estimate alpha"
#> [1] "Estimate q1"
#> [1] "Estimate q2"
#> [1] "Estimate q3"
#> [1] "Estimate beta"
#> [1] 9.768707e-10 1.974637e-09 0.000000e+00 1.820308e-09 9.465215e-03 1.596320e-07 1.088319e-03 3.839102e-04
#> [1] "Iteration 6"
#> [1] "Estimate p1"
#> [1] "Estimate p2"
#> [1] "Estimate p3"
#> [1] "Estimate alpha"
#> [1] "Estimate q1"
#> [1] "Estimate q2"
#> [1] "Estimate q3"
#> [1] "Estimate beta"
#> [1] 9.862053e-10 1.614438e-09 0.000000e+00 1.704341e-09 6.021229e-05 1.992563e-07 2.387598e-03 6.342837e-04
#> [1] "Iteration 7"
#> [1] "Estimate p1"
#> [1] "Estimate p2"
#> [1] "Estimate p3"
#> [1] "Estimate alpha"
#> [1] "Estimate q1"
#> [1] "Estimate q2"
#> [1] "Estimate q3"
#> [1] "Estimate beta"
#> [1]  5.177474e-08  1.874300e-09 -1.796250e+04  4.138904e-10  7.802773e-05  1.509075e-05  1.848810e-05  7.360461e-04
#> [1] "Iteration 8"
#> [1] "Estimate p1"
#> [1] "Estimate p2"
#> [1] "Estimate p3"
#> [1] "Estimate alpha"
#> [1] "Estimate q1"
#> [1] "Estimate q2"
#> [1] "Estimate q3"
#> [1] "Estimate beta"
#> [1]  1.059174e-08  1.484049e-07 -3.790187e+05  5.648965e-09  6.346873e-05  1.285139e-05  2.541129e-05  1.091563e-04
#> [1] "Iteration 9"
#> [1] "Estimate p1"
#> [1] "Estimate p2"
#> [1] "Estimate p3"
#> [1] "Estimate alpha"
#> [1] "Estimate q1"
#> [1] "Estimate q2"
#> [1] "Estimate q3"
#> [1] "Estimate beta"
#> [1]  6.296022e-06  3.229316e-08 -7.983213e+04  6.415143e-06  6.211295e-05  1.159249e-06  5.646427e-05  5.845998e-05

y_reconstructed <- array(0, c(nrow(p1), nrow(p2), nrow(p3), length(t_horizon)))
for (t in seq_len(length(t_horizon))) {
  for (j in seq_len(r)) {
    y_reconstructed[, , , t] <- y_reconstructed[, , , t] + hhat[[j]](t_horizon[t], j, kappa, result$alpha, bb, nu, mm) * result$p1[, j] %o% result$p2[, j] %o% result$p3[, j]
  }
  y_reconstructed[, , , t] <- y_reconstructed[, , , t] + ghat[[s4[t]]](t, s4[t], kappa, result$beta, aa, omega, mm) * result$q1[s1] %o% result$q2[s2] %o% result$q3[s3]
}




























































# Real data
library(doRNG)
load("data30buy.RData")
set.seed(1)

# 157203 x 672
print(dim(Process_buy))
fns <- c("length", "sum", "mean", "sd")
sample_count <- 100
sample_index_100 <- sample.int(nrow(Process_buy), size = sample_count, replace = FALSE)

sampled_data <- Process_buy[sample_index_100, 1:200]

print(dim(sampled_data))

all_times <- NULL
for (i in seq_len(nrow(sampled_data))) {
  for (j in seq_len(ncol(sampled_data))) {
    if (length(sampled_data[i, j][[1]]) > 0) {
      all_times <- c(all_times, sampled_data[i, j][[1]])
    }
  }
}

length_data <- matrix(NA, nrow(sampled_data), ncol(sampled_data))
for (i in seq_len(nrow(sampled_data))) {
  for (j in seq_len(ncol(sampled_data))) {
    length_data[i, j] <- length(x = sampled_data[i, j][[1]])
  }
}

all_times <- sort(unique(all_times))
print(length(all_times))

d <- 2
n1 <- 100
n2 <- 200
m1 <- 3
m2 <- 1
m3 <- 3
r <- 3
s1 <- kmeans(length_data, m1)$cluster
s2 <- rep(1, n2)
s3 <- c(rep(1, sum(all_times < 1/3)), rep(2, sum(all_times >= 1/3 & all_times < 2/3)), rep(3, sum(all_times >= 2/3)))
# t_horizon <- all_times
t_horizon <- 0:99 / 100
tt <- length(t_horizon)

y <- array(0, c(n1, n2, tt))
for (i in seq_len(nrow(sampled_data))) {
  for (j in seq_len(ncol(sampled_data))) {
    if (length(sampled_data[i, j][[1]]) > 0) {
      y[i, j, floor(sampled_data[i, j][[1]] * 100) + 1] <- unname(table(sampled_data[i, j][[1]]))
    }
  }
}

an <- 2
kappa <- 2
nu <- matrix(c(rep(mean(all_times < 1/3), r), rep(mean(all_times < 2/3), r)), r, an)
omega <- matrix(c(rep(mean(all_times < 1/3), m3), rep(mean(all_times < 2/3), m3)), m3, an)

mm <- an + kappa + 1

bb <- rep(list(rep(list(NA), mm)), r)
aa <- rep(list(rep(list(NA), mm)), m3)

for (j in seq_len(r)) {
  for (i in seq_len(mm)) {
    bb[[j]][[i]] <- function(t, j, i, kappa, nu) {
      if (i <= kappa + 1) {
        t^(i - 1)
      } else {
        max(t - nu[j, i - kappa - 1], 0)^kappa
      }
    }
  }
}

for (j in seq_len(m3)) {
  for (i in seq_len(mm)) {
    aa[[j]][[i]] <- function(t, j, i, kappa, omega) {
      if (i <= kappa + 1) {
        t^(i - 1)
      } else {
        max(t - omega[j, i - kappa - 1], 0)^kappa
      }
    }
  }
}

p1 <- matrix(rnorm(n1 * r), n1, r)
p2 <- matrix(rnorm(n2 * r), n2, r)
q1 <- seq_len(m1)
q2 <- seq_len(m2)
alpha <- matrix(0, r, mm)
beta <- matrix(0, m3, mm)

hhat <- rep(list(NA), r)
for (j in seq_len(r)) {
  hhat[[j]] <- function(t, j, kappa, alpha, bb, nu, mm) {
    sum(sapply(seq_len(mm), function(i) {
      alpha[j, i] * bb[[j]][[i]](t, j, i, kappa, nu)
    }))
  }
}

ghat <- rep(list(NA), m3)
for (j in seq_len(m3)) {
  ghat[[j]] <- function(t, j, kappa, beta, aa, omega, mm) {
    sum(sapply(seq_len(mm), function(i) {
      beta[j, i] * aa[[j]][[i]](t, j, i, kappa, omega)
    }))
  }
}

likelihood <- function(y, r, kappa, lambda, p1, p2, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) {
  yhat <- array(0, c(nrow(p1), nrow(p2), tt))
  for (t in seq_len(tt)) {
    for (j in seq_len(r)) {
      yhat[, , t] <- yhat[, , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j]
    }
    yhat[, , t] <- yhat[, , t] + ghat[[s3[t]]](t, s3[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2]
  }
  indicator <- y != 0
  omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[3]), function(t) indicator[, , t]))
  ll <- 0
  for (i1 in seq_len(nrow(p1))) {
    for (i2 in seq_len(nrow(p2))) {
      if (omega_set[i1, i2]) {
        ll <- ll + (y[i1, i2, ] - yhat[i1, i2, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, ] - yhat[i1, i2, ])
      }
    }
  }
  ll <- ll + lambda * (norm(p1, "F")^2 + norm(p2, "F")^2 + norm(q1[s1], "2")^2 + norm(q2[s2], "2")^2 + norm(alpha, "F")^2 + norm(beta, "F")^2)
  ll
}

dtrs <- function(y, r, kappa, lambda, p1, p2, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, epsilon = 1e-4) {
  iteration_l <- 0
  error <- 1
  while (error >= epsilon) {
    iteration_l <- iteration_l + 1
    print(paste("Iteration", iteration_l))
    # Estimate p1
    print("Estimate p1")
    cl <- parallel::makeCluster(parallel::detectCores(), outfile = "log/p1.log")
    doParallel::registerDoParallel(cl)
    i1 <- NULL
    p1new <- doRNG::`%dorng%`(foreach::foreach(i1 = seq_len(nrow(p1)), .combine = "rbind"), {
      optim(
        par = p1[i1, ],
        fn = function(p1, p2, q1, q2, s1, s2, s3, t_horizon, r, y, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
          yhat <- array(0, c(nrow(p2), length(t_horizon)))
          for (t in seq_len(length(t_horizon))) {
            for (j in seq_len(r)) {
              yhat[, t] <- yhat[, t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[j] * p2[, j]
            }
            yhat[, t] <- yhat[, t] + ghat[[s3[t]]](t, s3[t], kappa, beta, aa, omega, mm) * q1[s1][i1] * q2[s2]
          }
          indicator <- y != 0
          omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[2]), function(x) indicator[ , x]))
          result <- 0
          for (i2 in seq_len(nrow(p2))) {
            if (omega_set[i2]) {
              result <- result + (y[i2, ] - yhat[i2, ]) %*% diag(length(t_horizon)) %*% (y[i2, ] - yhat[i2, ])
            }
          }
          result <- result + lambda * norm(p1, type = "2")^2
        },
        p2 = p2,
        q1 = q1,
        q2 = q2,
        s1 = s1,
        s2 = s2,
        s3 = s3,
        t_horizon = t_horizon,
        r = r,
        y = y[i1, , ],
        alpha = alpha,
        beta = beta,
        hhat = hhat,
        ghat = ghat,
        aa = aa,
        bb = bb,
        mm = mm,
        nu = nu,
        omega = omega,
        lambda = lambda
      )$par
    })
    parallel::stopCluster(cl)
    jp1 <- 1 - likelihood(y, r, kappa, lambda, p1new, p2, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)
    # Estimate p2
    print("Estimate p2")
    cl <- parallel::makeCluster(parallel::detectCores(), outfile = "log/p2.log")
    doParallel::registerDoParallel(cl)
    i2 <- NULL
    p2new <- doRNG::`%dorng%`(foreach::foreach(i2 = seq_len(nrow(p2)), .combine = "rbind"), {
      optim(
        par = p2[i2, ],
        fn = function(p2, p1, q1, q2, s1, s2, s3, t_horizon, r, y, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
          yhat <- array(0, c(nrow(p1), length(t_horizon)))
          for (t in seq_len(length(t_horizon))) {
            for (j in seq_len(r)) {
              yhat[, t] <- yhat[, t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p2[j] * p1[, j]
            }
            yhat[, t] <- yhat[, t] + ghat[[s3[t]]](t, s3[t], kappa, beta, aa, omega, mm) * q2[s2][i2] * q1[s1]
          }
          indicator <- y != 0
          omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[2]), function(x) indicator[ , x]))
          result <- 0
          for (i1 in seq_len(nrow(p1))) {
            if (omega_set[i1]) {
              result <- result + (y[i1, ] - yhat[i1, ]) %*% diag(length(t_horizon)) %*% (y[i1, ] - yhat[i1, ])
            }
          }
          result <- result + lambda * norm(p2, type = "2")^2
        },
        p1 = p1,
        q1 = q1,
        q2 = q2,
        s1 = s1,
        s2 = s2,
        s3 = s3,
        t_horizon = t_horizon,
        r = r,
        y = y[, i2, ],
        alpha = alpha,
        beta = beta,
        hhat = hhat,
        ghat = ghat,
        aa = aa,
        bb = bb,
        mm = mm,
        nu = nu,
        omega = omega,
        lambda = lambda
      )$par
    })
    parallel::stopCluster(cl)
    jp2 <- 1 - likelihood(y, r, kappa, lambda, p1, p2new, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)
    # Estimate alpha
    print("Estimate alpha")
    alphanew <- optim(
      par = c(alpha),
      fn = function(alpha, p1, p2, q1, q2, s1, s2, s3, t_horizon, r, y, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
        alpha <- matrix(alpha, ncol = mm)
        yhat <- array(0, c(nrow(p1), nrow(p2), length(t_horizon)))
        for (t in seq_len(length(t_horizon))) {
          for (j in seq_len(r)) {
            yhat[, , t] <- yhat[, , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j]
          }
          yhat[, , t] <- yhat[, , t] + ghat[[s3[t]]](t, s3[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2]
        }
        indicator <- y != 0
        omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[3]), function(x) indicator[ , , x]))
        result <- 0
        for (i1 in seq_len(nrow(p1))) {
          for (i2 in seq_len(nrow(p2))) {
            if (omega_set[i1, i2]) {
              result <- result + (y[i1, i2, ] - yhat[i1, i2, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, ] - yhat[i1, i2, ])
            }
          }
        }
        result <- result + lambda * norm(alpha, type = "F")^2
      },
      p1 = p1,
      p2 = p2,
      q1 = q1,
      q2 = q2,
      s1 = s1,
      s2 = s2,
      s3 = s3,
      t_horizon = t_horizon,
      r = r,
      y = y,
      beta = beta,
      hhat = hhat,
      ghat = ghat,
      aa = aa,
      bb = bb,
      mm = mm,
      nu = nu,
      omega = omega,
      lambda = lambda
    )$par
    alphanew <- matrix(alphanew, ncol = mm)
    jalpha <- 1 - likelihood(y, r, kappa, lambda, p1, p2, q1, q2, s1, s2, s3, t_horizon, alphanew, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)

    jmax_index <- which.max(c(jp1, jp2, jalpha))
    if (jmax_index == 1) {
      p1 <- p1new
    } else if (jmax_index == 2) {
      p2 <- p2new
    } else if (jmax_index == 3) {
      alpha <- alphanew
    }
    # Estimate q1
    print("Estimate q1")
    q1new <- optim(
      par = q1,
      fn = function(q1, p1, p2, q2, s1, s2, s3, t_horizon, r, y, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
        yhat <- array(0, c(nrow(p1), nrow(p2), length(t_horizon)))
        for (t in seq_len(length(t_horizon))) {
          for (j in seq_len(r)) {
            yhat[, , t] <- yhat[, , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j]
          }
          yhat[, , t] <- yhat[, , t] + ghat[[s3[t]]](t, s3[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2]
        }
        indicator <- y != 0
        omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[3]), function(x) indicator[ , , x]))
        result <- 0
        for (i1 in seq_len(nrow(p1))) {
          for (i2 in seq_len(nrow(p2))) {
            if (omega_set[i1, i2]) {
              result <- result + (y[i1, i2, ] - yhat[i1, i2, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, ] - yhat[i1, i2, ])
            }
          }
        }
        result <- result + lambda * norm(q1, type = "2")^2
      },
      p1 = p1,
      p2 = p2,
      q2 = q2,
      s1 = s1,
      s2 = s2,
      s3 = s3,
      t_horizon = t_horizon,
      r = r,
      y = y,
      alpha = alpha,
      beta = beta,
      hhat = hhat,
      ghat = ghat,
      aa = aa,
      bb = bb,
      mm = mm,
      nu = nu,
      omega = omega,
      lambda = lambda
    )$par
    jq1 <- 1 - likelihood(y, r, kappa, lambda, p1, p2, q1new, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)
    # Estimate q2
    print("Estimate q2")
    q2new <- optim(
      par = q2,
      fn = function(q2, p1, p2, q1, s1, s2, s3, t_horizon, r, y, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
        yhat <- array(0, c(nrow(p1), nrow(p2), length(t_horizon)))
        for (t in seq_len(length(t_horizon))) {
          for (j in seq_len(r)) {
            yhat[, , t] <- yhat[, , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j]
          }
          yhat[, , t] <- yhat[, , t] + ghat[[s3[t]]](t, s3[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2]
        }
        indicator <- y != 0
        omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[3]), function(x) indicator[ , , x]))
        result <- 0
        for (i1 in seq_len(nrow(p1))) {
          for (i2 in seq_len(nrow(p2))) {
            if (omega_set[i1, i2]) {
              result <- result + (y[i1, i2, ] - yhat[i1, i2, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, ] - yhat[i1, i2, ])
            }
          }
        }
        result <- result + lambda * norm(q2, type = "2")^2
      },
      p1 = p1,
      p2 = p2,
      q1 = q1,
      s1 = s1,
      s2 = s2,
      s3 = s3,
      t_horizon = t_horizon,
      r = r,
      y = y,
      alpha = alpha,
      beta = beta,
      hhat = hhat,
      ghat = ghat,
      aa = aa,
      bb = bb,
      mm = mm,
      nu = nu,
      omega = omega,
      lambda = lambda
    )$par
    jq2 <- 1 - likelihood(y, r, kappa, lambda, p1, p2, q1, q2new, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)
    # Estimate beta
    print("Estimate beta")
    betanew <- optim(
      par = c(beta),
      fn = function(beta, p1, p2, q1, q2, s1, s2, s3, t_horizon, r, y, alpha, hhat, ghat, aa, bb, mm, nu, omega, lambda) {
        beta <- matrix(beta, ncol = mm)
        yhat <- array(0, c(nrow(p1), nrow(p2), length(t_horizon)))
        for (t in seq_len(length(t_horizon))) {
          for (j in seq_len(r)) {
            yhat[, , t] <- yhat[, , t] + hhat[[j]](t_horizon[t], j, kappa, alpha, bb, nu, mm) * p1[, j] %o% p2[, j]
          }
          yhat[, , t] <- yhat[, , t] + ghat[[s3[t]]](t, s3[t], kappa, beta, aa, omega, mm) * q1[s1] %o% q2[s2]
        }
        indicator <- y != 0
        omega_set <- Reduce(`|`, lapply(seq(dim(indicator)[3]), function(x) indicator[ , , x]))
        result <- 0
        for (i1 in seq_len(nrow(p1))) {
          for (i2 in seq_len(nrow(p2))) {
            if (omega_set[i1, i2]) {
              result <- result + (y[i1, i2, ] - yhat[i1, i2, ]) %*% diag(length(t_horizon)) %*% (y[i1, i2, ] - yhat[i1, i2, ])
            }
          }
        }
        result <- result + lambda * norm(beta, type = "2")^2
      },
      p1 = p1,
      p2 = p2,
      q1 = q1,
      q2 = q2,
      s1 = s1,
      s2 = s2,
      s3 = s3,
      t_horizon = t_horizon,
      r = r,
      y = y,
      alpha = alpha,
      hhat = hhat,
      ghat = ghat,
      aa = aa,
      bb = bb,
      mm = mm,
      nu = nu,
      omega = omega,
      lambda = lambda
    )$par
    betanew <- matrix(betanew, ncol = mm)
    jbeta <- 1 - likelihood(y, r, kappa, lambda, p1, p2, q1, q2, s1, s2, s3, t_horizon, alpha, betanew, hhat, ghat, aa, bb, mm, nu, omega) / likelihood(y, r, kappa, lambda, p1, p2, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega)

    jmax_index <- which.max(c(jq1, jq2, jbeta))
    if (jmax_index == 1) {
      q1 <- q1new
    } else if (jmax_index == 2) {
      q2 <- q2new
    } else if (jmax_index == 3) {
      beta <- betanew
    }
    error <- max(jp1, jp2, jalpha, jq1, jq2, jbeta)
    print(c(jp1, jp2, jalpha, jq1, jq2, jbeta))
    saveRDS(list(p1 = p1, p2 = p2, q1 = q1, q2 = q2, alpha = alpha, beta = beta), paste0("dtrsreal", iteration_l, ".rds"))
  }
  return(list(p1 = p1, p2 = p2, q1 = q1, q2 = q2, alpha = alpha, beta = beta))
}

result <- dtrs(y, r, kappa, lambda = 5, p1, p2, q1, q2, s1, s2, s3, t_horizon, alpha, beta, hhat, ghat, aa, bb, mm, nu, omega, epsilon = 1e-4)
# result <- readRDS("dtrsreal3.rds")

# Plot the results
x <- seq(0, 1, length.out = 100)
molten <- data.frame(x = NULL, y = NULL, index = NULL, label = NULL)
for (j in seq_len(r)) {
  for (data_i in seq_len(length(x))) {
    molten <- rbind(molten, data.frame(x = x[data_i], y = hhat[[j]](x[data_i], j, kappa, result$alpha, bb, nu, mm), index = j, label = "h"))
  }
}
for (j in seq_len(m3)) {
  for (data_i in seq_len(length(x))) {
    molten <- rbind(molten, data.frame(x = x[data_i], y = ghat[[j]](x[data_i], j, kappa, result$beta, aa, omega, mm), index = j, label = "g"))
  }
}

pdf("dtrsreal.pdf", width = 8, height = 6)
ggplot2::ggplot(molten, ggplot2::aes(x = x, y = y, color = factor(index))) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(~ label, ncol = 1) +
  ggplot2::theme_bw()
dev.off()

y_reconstructed <- array(0, c(nrow(p1), nrow(p2), length(t_horizon)))
for (t in seq_len(length(t_horizon))) {
  for (j in seq_len(r)) {
    y_reconstructed[, , t] <- y_reconstructed[, , t] + hhat[[j]](t_horizon[t], j, kappa, result$alpha, bb, nu, mm) * result$p1[, j] %o% result$p2[, j]
  }
  y_reconstructed[, , t] <- y_reconstructed[, , t] + ghat[[s3[t]]](t, s3[t], kappa, result$beta, aa, omega, mm) * result$q1[s1] %o% result$q2[s2]
}

sampled_data <- Process_buy[sample_index_100, ]
all_trading_activities_per_day <- rep(list(NULL), ncol(sampled_data))
for (j in seq_len(ncol(sampled_data))) {
  for (i in seq_len(nrow(sampled_data))) {
    if (length(sampled_data[i, j][[1]]) > 0) {
      all_trading_activities_per_day[[j]] <- c(all_trading_activities_per_day[[j]], sampled_data[i, j][[1]])
    }
  }
}

sampled_mean_data_transition <- list()
for (day in seq_len(ncol(sampled_data))) {
  day_density <- density(all_trading_activities_per_day[[day]], na.rm = TRUE)
  sampled_mean_data_transition[[day]] <- list(
    visible = FALSE,
    name = paste0("day = ", day),
    x = day_density$x,
    y = day_density$y
  )
}
sampled_mean_data_transition[3][[1]]$visible <- TRUE

steps <- list()
fig <- plotly::plot_ly()
for (i in seq_len(ncol(sampled_data))) {
  fig <- plotly::add_lines(fig,
    x = sampled_mean_data_transition[i][[1]]$x, y = sampled_mean_data_transition[i][[1]]$y, visible = sampled_mean_data_transition[i][[1]]$visible,
    name = sampled_mean_data_transition[i][[1]]$name, type = "scatter", mode = "lines", hoverinfo = "name", fill = "tozeroy",
    line = list(color = "00CED1"), showlegend = FALSE
  )

  day <- list(
    args = list("visible", rep(FALSE, length(sampled_mean_data_transition))),
    method = "restyle"
  )
  day$args[[2]][i] <- TRUE
  steps[[i]] <- day
}

# add slider control to plot
fig <- plotly::layout(fig, sliders = list(list(
    active = 3,
    currentvalue = list(prefix = NULL),
    steps = steps
  )))

htmlwidgets::saveWidget(
  widget = fig, #the plotly object
  file = "3.html", #the path & file name
  selfcontained = TRUE #creates a single html file
)
