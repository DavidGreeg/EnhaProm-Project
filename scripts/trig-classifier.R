# library(mosaic)
library(rootSolve)
coss <- function(case) {
  return(switch(case,
    expression(cos(1 * pi * x)),
    expression(cos(2 * pi * x)),
    expression(cos(3 * pi * x)),
    expression(cos(4 * pi * x)),
    expression(cos(5 * pi * x)),
    expression(cos(6 * pi * x)),
    expression(cos(7 * pi * x)),
    expression(cos(8 * pi * x)),
  ))
}

sins <- function(case) {
  return(switch(case,
    expression(sin(1 * pi * x)),
    expression(sin(2 * pi * x)),
    expression(sin(3 * pi * x)),
    expression(sin(4 * pi * x)),
    expression(sin(5 * pi * x)),
    expression(sin(6 * pi * x)),
    expression(sin(7 * pi * x)),
    expression(sin(8 * pi * x)),
  ))
}

cossins <- function(case) {
  return(switch(case,
    expression(cos(1 * pi * x)),
    expression(sin(1 * pi * x)),
    expression(cos(2 * pi * x)),
    expression(sin(2 * pi * x)),
    expression(cos(3 * pi * x)),
    expression(sin(3 * pi * x)),
    expression(cos(4 * pi * x)),
    expression(sin(4 * pi * x)),
    expression(cos(5 * pi * x)),
    expression(sin(5 * pi * x)),
    expression(cos(6 * pi * x)),
    expression(sin(6 * pi * x)),
    expression(cos(7 * pi * x)),
    expression(sin(7 * pi * x)),
    expression(cos(8 * pi * x)),
    expression(sin(8 * pi * x)),
  ))
}

comb_cases <- function(k = 4, coss = FALSE, sins = FALSE,
                       print = TRUE) {
  if (coss || sins) {
    set <- 1:k
    cases <- c()
    peaks <- c()
    case_names <- c()
    for (i in set) {
      combs <- combn(set, i)
      combs_len <- length(combs[1, ])
      for (j in 1:combs_len) {
        if (coss && sins) {
          combcase <- cases_set(combs[, j], coss = TRUE, sins = TRUE)
          case_name <- case_names_set(combs[, j], coss = TRUE, sins = TRUE)
        } else {
          if (coss) {
            combcase <- cases_set(combs[, j], coss = TRUE)
            case_name <- case_names_set(combs[, j], coss = TRUE)
          }
          if (sins) {
            combcase <- cases_set(combs[, j], sins = TRUE)
            case_name <- case_names_set(combs[, j], sins = TRUE)
          }
        }
        case_stri <- as.character(parse(text = combcase))
        eva_fun <- as.function(alist(x = , eval(D(combcase, "x"))))
        num_peaks <- length(uniroot.all(eva_fun, c(0, 1)))
        if (print)
          cat("CASE: ", case_stri, "\t", num_peaks, "\n")
        cases <- c(cases, case_stri)
        peaks <- c(peaks, num_peaks)
        case_names <- c(case_names, case_name)
      }
    }
    return(data.frame(caseExpr = cases, peakNum = peaks, caseName = case_names))
  } else {
    stop("No function option selected")
  }
}

conc <- function(x, y) parse(text = paste(x, "+", y))
# : eval_func <- function(ex) function(x) eval(D(ex, "x"))

cases_set <- function(set, coss = FALSE, sins = FALSE) {
  if (coss || sins) {
    flag <- FALSE
    for (elem in set) {
      if (coss && sins) {
        caseset <- ifelse(flag, conc(caseset, cossins(elem)), cossins(elem))
      } else {
        if (coss)
          caseset <- ifelse(flag, conc(caseset, coss(elem)), coss(elem))
        if (sins)
          caseset <- ifelse(flag, conc(caseset, sins(elem)), sins(elem))
      }
      if (!flag)
        flag <- TRUE
    }
    return(caseset)
  } else {
    stop("No function option selected")
  }
}

combs_no_rep <- function(set, k) {
  return(combn(set, k, simplify = TRUE))
}

case_names_set <- function(set, coss = FALSE, sins = FALSE) {
  if (coss || sins) {
    flag <- FALSE
    for (elem in set) {
      if (coss && sins) {
        caseset <- ifelse(flag,
                          paste(caseset, cossins_names(elem), sep = "+"),
                          cossins_names(elem))
      } else {
        if (coss)
          caseset <- ifelse(flag,
                            paste(caseset, coss_names(elem), sep = "+"),
                            coss_names(elem))
        if (sins)
          caseset <- ifelse(flag,
                            paste(caseset, sins_names(elem), sep = "+"),
                            sins_names(elem))
      }
      if (!flag)
        flag <- TRUE
    }
    return(caseset)
  } else {
    stop("No function option selected")
  }
}

coss_names <- function(case) {
  return(switch(case,
    "c1", "c2", "c3", "c4",
    "c5", "c6", "c7", "c8",
  ))
}

sins_names <- function(case) {
  return(switch(case,
    "s1", "s2", "s3", "s4",
    "s5", "s6", "s7", "s8",
  ))
}

cossins_names <- function(case) {
  return(switch(case,
    "c1", "s1", "c2", "s2", "c3", "s3", "c4", "s4",
    "c5", "s5", "c6", "s6", "c7", "s7", "c8", "s8",
  ))
}

# Function to create and evaluate sequential lm models
sequential_lm1 <- function(data, response, predictor,
                           max_terms, bottomtop = FALSE) {
  models <- list() # List to store the resulting models

  # To do later:
  # - function that gets the relative number of peaks by
  #   looking non-increasing values in the response variable.

  # Develop this mechanism later
  # NOTE: Try this using 'ifelse'
  if (!bottomtop) {
    seq_terms <- seq(max_terms, 1, by = -1)
  } else {
    seq_terms <- seq(1, max_terms, by = 1)
  }

  for (n in seq_terms) {
    if (n < max_terms) {
      # Add the current term to the formula
      formula_terms <- paste(formula_terms, "+",
                             paste0("cos(", n, " * 2 * pi * ", predictor, ")"))
    } else {
      # Start with the base formula
      formula_terms <- paste0("cos(", max_terms, " * 2 * pi * ", predictor, ")")
    }

    # Create the full formula
    formula <- as.formula(paste(response, "~", formula_terms))

    # Fit the model
    model <- lm(formula, data = data)

    # Store the model in the list
    models[[paste0("lm_", max_terms - n + 1)]] <- model
  }
  return(models)
}
# Example usage
# Assuming `my_data` is a data frame with columns `spc_len` (response) and `spc_pos` (predictor):
# models <- sequential_lm(my_data, response = "spc_len", predictor = "spc_pos", max_terms = 5)

compare_cases <- function(list_functions) return(list_functions)

sequential_lm <- function(data, response, predictor, n_peaks, silent = FALSE,
                          bottomtop = FALSE, aem_option = 1, plot = TRUE,
                          more_sliders = TRUE, even_more_sliders = FALSE,
                          check_again = FALSE, check_again_n = 0) {
  #: models <- list() # List to store the resulting models

  # To do later:
  # - function that gets the relative number of peaks by
  #   looking non-increasing values in the response variable
  #   (after removing ~80% of lower-end points).

  # Develop this mechanism later
  # NOTE: Try this using 'ifelse'
  if (!bottomtop) {
    i_peaks <- seq(n_peaks * 2, 2, by = -2)
    peakstart <- (n_peaks * 2)
  } else {
    i_peaks <- seq(2, n_peaks * 2, by = 2)
    peakstart <- 2
  }

  # Default model
  #: f_default <- paste0(response, " ~ cos(pi * ", predictor, ")")
  #: best_aem <- adj_error_metric(lm(as.formula(f_default), data = data))
  best_aem <- -22
  max_aem <- best_aem
  #: print(max_aem)

  all_done <- FALSE
  cntimes <- check_again_n
  if (check_again) cntimes <- 1
  if (cntimes) i_peaks <- c(i_peaks, rep(0, cntimes))
  cstorage_size <- ifelse(cntimes,
                          length(i_peaks) - cntimes, length(i_peaks))

  # Vectors where we'll store coefficients
  dim_arr <- 2 * cstorage_size
  coeffs_storage <- array(rep(0, 2 * dim_arr), dim = c(dim_arr, 2),
                          dimnames = list(1:dim_arr, c("l_coefs", "s_coefs")))
  #: slider_storage <- c()
  pcoef_tosave <- 0
  scoef_tosave <- 0

  slide_floats <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
  if (more_sliders)
    slide_floats <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
  if (even_more_sliders)
    slide_floats <- seq(0, 0.5, by = 0.025)
  slide_terms <- as.character(slide_floats[-1])
  slide_terms <- sapply(slide_terms, function(x) paste0(" - ", x, "))"))
  slide_terms <- c("))", slide_terms)
  for (ipeak in i_peaks) {
    if (!bottomtop) {
      pcoefs <- ipeak:(ipeak - 2)
      if (!(ipeak - 2)) pcoefs <- pcoefs[-length(pcoefs)]
      #  ^This removes the case of 'ipeak == 0'
    } else {
      pcoefs <- (ipeak - 2):ipeak
      #:::::::::::::::::: if (!ipeak) pcoefs <- pcoefs[-1]
      if (!silent) cat("ipeak: ", ipeak, "peakstart: ", peakstart, "\n")
    }

    # Only used if 'check_again=TRUE'
    # One last run to try non-used "pi_coefs"
    if (!ipeak) {
      if (all_done) print("No further checkings necessary")
      if (all_done) break
      pcoefs <- which(coeffs_storage[, "l_coefs"] == 0)
      # print("I'm in uwu")
      # print(pcoefs)
      # all_done <- TRUE
    }

    #: p_models <- list()
    i <- 0
    for (pcoef in pcoefs) {
      i <- i + 1
      if (pcoef == pcoef_tosave) next
      if (!ipeak) all_done <- TRUE
      if (ipeak != peakstart) {
        # f_terms <- paste(last_best_terms, "+",
        f_terms <- paste(last_f, "+",
                         paste0("cos(pi * ((", pcoef, " * ", predictor, ")))"))
      } else {
        f_terms <- paste0("cos(pi * ((", pcoef, " * ", predictor, ")))")
      }
      t_list <- sapply(slide_terms,
                       function(s_term) gsub("))$", s_term, f_terms))
      f_list <- paste(response, "~", t_list)
      lm_list <- lapply(f_list,
                        function(f_i) lm(as.formula(f_i), data = data))
      aem_list <- sapply(lm_list,
                         function(i_lm) adj_error_metric(i_lm))

      newmax_aem <- max(aem_list)

      f_index <- which(aem_list == newmax_aem)

      if (!silent) {
        cat("\n==============================|==============================\n")
        print(aem_list)
        print(f_index)
        print(slide_floats[f_index])
      }
      # print(lm_list[[f_index]])

      if (newmax_aem > max_aem) {
        if (!silent) cat("New max AEM!\n")
        if (!ipeak) all_done <- FALSE
        # ^This conditional ('all_done') will help us break the loop
        # if no new terms are added while "checking again"

        pcoef_tosave <- pcoef

        scoef_tosave <- slide_floats[f_index]

        newbest_lm <- lm_list[[f_index]]
        lm_coefficients <- newbest_lm$coefficients
        lcoef_tosave <- lm_coefficients[length(lm_coefficients)]

        max_aem <- newmax_aem

        selected_t <- t_list[f_index]
        #: ind_tostore <- abs(ipeak - pcoef) + 1
        #: inds_storage <- c(0, 0, 0)
        #: inds_storage[ind_tostore] <- 1
      }

      #: f_best <- compare_cases(f_list) ###
      #: f_terms <- f_best[1]
      #: f_model <- f_best[2]

      # Save 'f_best' (best model in loop)
      #: f_pcoef <- as.formula(paste(response, "~", f_terms))
      #: f_model <- lm(f_pcoef, data = data)
      #: p_models[[paste0("lm_", i)]] <- f_model
      if (!silent) {
        cat(t_list, sep = "\n")
        cat("Max AEM: ", max_aem, "; pcoef: ", pcoef_tosave,"; scoef: ", scoef_tosave, "\n")
      }
    }
    if (max_aem > best_aem) {
      best_aem <- max_aem
      best_lm <- newbest_lm
      last_f <- selected_t
      coeffs_storage[pcoef_tosave, "l_coefs"] <- lcoef_tosave
      coeffs_storage[pcoef_tosave, "s_coefs"] <- scoef_tosave
      if (!silent) {
        print(coeffs_storage[, "l_coefs"])
        print(coeffs_storage[, "s_coefs"])
      }
    }

    #: last_f <- selected_f
    #: models[[paste0("lms_", ipeak)]] <- p_models
  }
  coeffs_inds <- coeffs_storage[, "l_coefs"] != 0
  coeffs_storage[coeffs_inds, "l_coefs"] <- best_lm$coefficients[-1]

  # Visualize final array of pi and slide coefficients
  if (silent) {
    print(t(coeffs_storage))
  } else {
    print(coeffs_storage)
  }

  if (plot) {
    fitdata <- data.frame(predictor = seq(min(data[, predictor]),
                                          max(data[, predictor]),
                                          length.out = 250))
    colnames(fitdata)[1] <- predictor
    fitdata$pred_response <- predict(best_lm, newdata = fitdata)

    plot(y = data[, response], x = data[, predictor])
    points(data[, predictor], fitted(best_lm), col = "red", pch = 20)
    lines(x = fitdata[, predictor], y = fitdata$pred_response)
  }
  return(best_lm)
}

source("/home/davidfm/Projects/UBMI-IFC/EnhaProm/scripts/genome-functions.R")
lm_pvalue <- lm_pvalue

adj_error_metric <- function(lmodel, best_adj = FALSE,
                             silent = TRUE, brief = FALSE, aem_option = 1) {
  summ_lmodel <- summary(lmodel)

  adj_rsq <- summ_lmodel$adj.r.squared
  rsq <- summ_lmodel$r.squared
  res_err <- summ_lmodel$sigma
  p_val <- lm_pvalue(summ_lmodel)
  #: t_mean <- tvalues(summ_lmodel, mean = TRUE)
  #: t_prod <- tvalues(summ_lmodel, prod = TRUE)

  if (adj_rsq > 0.999 && rsq > 0.999) {
    if (!silent)
      cat("\t ~~>AEM:  OVER-FITTED")
    return(0.00001)
  }

  errs_names <- c("Adj. R-Squared", "Res. Std. Error",
                  "R-Squared", "P-Value")
  errs_vectr <- c(adj_rsq, res_err, rsq, p_val)
  #: errs_names <- c("Adj. R-Squared", "Res. Std. Error",
  #:                 "R-Squared", "P-Value", "Mean T", "Prod T")
  #: errs_vectr <- c(adj_rsq, res_err, rsq, p_val, t_mean, t_prod)

  aem <- switch(aem_option,
    (adj_rsq * res_err) / (rsq * p_val),
    (adj_rsq * p_val) / (rsq * res_err),
    adj_rsq / (rsq * p_val * res_err),
    res_err^2 / p_val^2,
    adj_rsq,
    rsq
  )
  names(errs_vectr) <- errs_names

  if (brief) {
    cat("\t ~~>AEM: ", aem)
    silent <- TRUE
  }
  if (!silent) {
    print(lmodel$coefficients)
    print("---  ----  ----  ----  ---  -··· ···-  ----  ----  ----  ----  ---")
    print(errs_vectr)
    cat("\n                     ~~>AEM: ", aem, "\n")
    print("=============================|||-|||==============================")
  }
  if (best_adj)
    return(adj_rsq)
  return(aem)
}



