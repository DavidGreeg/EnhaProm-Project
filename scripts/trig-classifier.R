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
