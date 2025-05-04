library(knitr)
library(kableExtra)

library(stringi)
library(stringr)

library(dplyr)

library(ggplot2)
library(cowplot)
library(see)

library(irlba)

#______________________________________________________________________________#
#==================================LIBRARIE====================================#
required_libs <- function(libs) {
  for (lib in libs)
  suppressPackageStartupMessages(library(lib, character.only = TRUE))
}

required_libs(c("knitr", "kableExtra", "stringi", "dplyr"))
# teatable requires those libs

#______________________________________________________________________________#
#===================================TEATABLE===================================#
teatable <- function(tabl, cat = FALSE, decay = TRUE, pale = TRUE, colsize) {

  if (!decay) {
    rowcolors <-  c("#a3cfa3", "#b3e6b3", "green!20")
  } else {
    if (!pale) {
      rowcolors <- c("#8c8c66", "#e6e6cc", "#c2c299")
    } else {
      rowcolors <- c("#ccccb3", "#f4f4e2", "#e0e0cc")
    }
  }

  kbl_tabl <- kbl(tabl, align = "c", booktabs = TRUE, longtable = FALSE) %>%
    kable_styling(position = "center",
                  latex_options = c("striped", "scale_down", "hold_position"),
                  stripe_color = rowcolors[3], full_width = FALSE) %>%
    row_spec(seq(2, nrow(tabl), by = 2), background = rowcolors[2]) %>%
    row_spec(0, background = rowcolors[1]) #%>%

  if (!missing(colsize)) {
    n_cols <- ncol(tabl)
    kbl_tabl <- kbl_tabl %>% column_spec(1:n_cols, width = colsize)
  }

  kbl_tabl <- gsub(x = kbl_tabl,
                   "\\\\midrule",
                   "\\\\specialrule\\{0.6pt\\}\\{0.8pt\\}\\{0.6pt\\}")
  kbl_tabl <- gsub(x = kbl_tabl,
                   "\\\\toprule",
                   "\\\\specialrule\\{1pt\\}\\{0pt\\}\\{1pt\\}")
  kbl_tabl <- gsub(x = kbl_tabl,
                   "\\\\bottomrule",
                   "\\\\specialrule\\{1pt\\}\\{0.6pt\\}\\{0pt\\}")
  kbl_tabl <- gsub(x = kbl_tabl,
                   "\\\\begin\\{table\\}([^\n]*)",
                   "\\\\begin\\{table\\}\\1\n\\\\footnotesize")

  if (!cat) {asis_output(kbl_tabl)} else {cat(kbl_tabl)}
}

rm_lastchar <- function(chain) return(substr(chain, 1, nchar(chain) - 1))

ncol_latex <- function(btab_str) {
  return(nchar(strsplit(rm_lastchar(btab_str), "\\}.*\\{")[[1]][2]))
  # nchar(gsub("\\{[a-z0-9.]+\\}|\\|", "",
  #            str_split_fixed(substr(btab_str, 1,
  #                                   nchar(btab_str)-1), "\\{", n=3)[3])))
  # ^more extensive/broader take on possible '\begin{tabular}' configurations
}

#______________________________________________________________________________#
#===============================STATE-MODIFIERS================================#
state_mod <- function(logvect, inds) {
  inds <- ifelse(inds & length(logvect) >= inds,
                 inds, stop("inds out of bounds or is NULL"))
  logvect[inds] <- !logvect[inds]
  return(logvect)
}

int_state_mod <- function(intvect, addvect, logvect, inds) {
  inds <- ifelse(inds & length(logvect) >= inds,
                 inds, stop("inds out of bounds or is NULL"))
  addvect <- ((-1)^logvect[inds]) * addvect
  intvect[inds] <- intvect[inds] + addvect
  return(intvect)
}

#______________________________________________________________________________#
#================================GET-FROM-TABL=================================#
get_rowspans <- function(lines_mrow) {
  return(as.numeric(gsub(x = lines_mrow,
                         ".*\\\\multirow\\{-?([0-9]+)\\}.*", "\\1")))
}

get_colsigns <- function(mrow_vec) {
  return(mrow_vec / replace(abs(mrow_vec), mrow_vec == 0, 1))
}

#______________________________________________________________________________#
#===============================COLOR-MULTIROWS================================#
color_multirows <- function(lines, ncols) {
  logic_vec <- rep(0, ncols)
  integ_vec <- rep(0, ncols)

  if (!length(lines)) {
    stop("No argument 'lines' provided")
  } else {
    len_lines <- length(lines)
  }

  search_str <- "\\[.*\\]\\{(.*) (.*)\\}"
  replace_str1 <- "\\{\\1\\\\cellcolor\\{brown!10\\} \\2\\}"
  replace_str2 <- "\\{\\1\\\\cellcolor\\{brown!25\\} \\2\\}"
  topaste_str1 <- "\\cellcolor{brown!10}"
  topaste_str2 <- "\\cellcolor{brown!25}"

  return(sapply(lines[len_lines:1], function(line) {
    line_cells <- unlist(strsplit(line, " & "))
    multirow_col_inds <- grep("\\\\multirow\\{", line_cells)

    if (length(multirow_col_inds)) {
      logic_vec <<- state_mod(logic_vec, multirow_col_inds)
      n_rowspans <- get_rowspans(line_cells[multirow_col_inds]) - 1
      integ_vec <<- int_state_mod(integ_vec, n_rowspans,
                                  logic_vec, multirow_col_inds)
      for (col_i in multirow_col_inds)
      line_cells[col_i] <- ifelse(logic_vec[col_i] == 0,
                                  gsub(search_str, replace_str1,
                                       line_cells[col_i]),
                                  gsub(search_str, replace_str2,
                                       line_cells[col_i]))
    } else if (sum(abs(integ_vec))) {
      for (col_i in (1:ncols)[as.logical(integ_vec)])
      line_cells[col_i] <- ifelse(logic_vec[col_i] == 0,
                                  paste0(topaste_str1, line_cells[col_i]),
                                  paste0(topaste_str2, line_cells[col_i]))
      integ_vec <<- integ_vec - get_colsigns(integ_vec)
    }

    paste(line_cells, collapse = " & ")
  })[len_lines:1])
}

#______________________________________________________________________________#
#=================================COFFEETABLE==================================#
coffeetable <- function(tabl, collapse = TRUE, cat = FALSE,
                        row.names = TRUE, row.height = 1.3) {
  # The name is actually a reference to the amount of coffee needed
  #   to get to this workaround solution, 5 seconds later I noticed
  #   the ironic coincidence between 'tea' and 'coffee'.
  #   Note: this was before I changed the set to a brown-scale
  if (!row.names)
    rownames(tabl) <- NULL

  kbl_tabl <- kbl(tabl, booktabs = TRUE, longtable = FALSE, linesep = "") %>%
    kable_styling(position = "center",
                  latex_options = c("scale_down", "hold_position")) #%>%

  if (collapse)
    kbl_tabl <- kbl_tabl %>% collapse_rows()

  tabl_str <- toString(kbl_tabl)
  tabl_str_lines <- unlist(strsplit(tabl_str, "\n"))

  if (collapse) {
    cmid_lines <- grep("\\\\cmidrule\\{", tabl_str_lines)
    tabl_str_lines <- tabl_str_lines[-cmid_lines]
  }

  rowcolor1 <- "\\rowcolor{brown!10}\n"
  rowcolor2 <- "\\rowcolor{brown!25}\n"
  rowcolor3 <- "\\rowcolor{brown!42}\n"

  rheight_str1 <- "\\renewcommand{\\arraystretch}{"
  rheight_str2 <- "\n\\renewcommand{\\arraystretch}{1}\n"

  tabu_lines <- grep("\\\\(?:begin|end)\\{tabular\\}", tabl_str_lines)
  tabl_str_lines[tabu_lines[1] - 1] <- paste0(rheight_str1,
                                              as.character(row.height), "}\n",
                                              tabl_str_lines[tabu_lines[1] - 1])
  tabl_str_lines[tabu_lines[1] - 1] <- paste0("\\footnotesize\n",
                                              tabl_str_lines[tabu_lines[1] - 1])
  tabl_str_lines[tabu_lines[2]] <- paste0(tabl_str_lines[tabu_lines[2]],
                                          rheight_str2)

  inter_lines <- grep("\\b\\\\midrule\\b|\\b\\\\bottomrule\\b", tabl_str_lines)
  all_i <- seq(inter_lines[1] + 1, inter_lines[2] - 1)
  evn_i <- seq(inter_lines[1] + 1, inter_lines[2] - 1, 2)
  odd_i <- all_i[! all_i %in% evn_i]

  if (collapse) {
    tabl_ncols <- ncol_latex(tabl_str_lines[tabu_lines[1]])
    tabl_str_lines[all_i] <- color_multirows(tabl_str_lines[all_i], tabl_ncols)
  }

  tabl_str_lines[all_i[1] - 2] <- paste0(rowcolor3, tabl_str_lines[all_i[1]-2])
  # Replace '\toprule' and '\midrule'.
  tabl_str_lines[all_i[1] - 3] <- "\\specialrule{1pt}{0pt}{0.5pt}"
  tabl_str_lines[all_i[1] - 1] <- "\\specialrule{0.6pt}{1pt}{0pt}"
  # This may also be possible by setting '\aboverulesep' and '\belowrulesep'

  tabl_str_lines[all_i] <- paste0(tabl_str_lines[all_i], "[-0.12pt]")
  tabl_str_lines[evn_i] <- paste0(rowcolor2, tabl_str_lines[evn_i])
  tabl_str_lines[odd_i] <- paste0(rowcolor1, tabl_str_lines[odd_i])
  tabl_str <- paste(tabl_str_lines, collapse = "\n")

  if (!cat) {asis_output(tabl_str)} else {cat(tabl_str)}
}

#______________________________________________________________________________#
#==============================GET-LEGEND-BYPASS===============================#
# NOTE: 'cowplot::get_legend' doesn't seem to work right with 'ggplot2=3.5*'.
# This solution was taken from: 'https://github.com/wilkelab/cowplot/issues/202'
get_legend_bypass <- function(plot) {
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  nonzero <- vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE)
  idx <- which(nonzero)
  if (length(idx) > 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}

#______________________________________________________________________________#
#==================================PLOTTING====================================#
Type <- factor()	# This variables are just 'decoys' so my object-linter
Means <- integer()	# doesn't fire a lot of warnings for the following
StDevs <- integer() # plotting functions
Field <- factor()

# To plot pyramid-plots
pyrplot_ <- function(cre_data, kmer_labels, title,
                     x_label, y_label, y_breaks = NULL) {
  CREs <- c("Enhancer", "Promoter")
  field_order <- filter(cre_data, Type == CREs[1])$Field
  fact_field_order <- factor(cre_data$Field, field_order)

  if (missing(kmer_labels)) kmer_labels <- field_order

  ggplot(cre_data) +
    geom_bar(aes(x = fact_field_order,
                 y = ifelse(Type == CREs[1], -Means, Means),
                 fill = paste(Type, "Means")),
             stat = "identity", position = "identity",
             alpha = 0.6, width = 0.7) +
    geom_errorbar(aes(x = fact_field_order,
                      ymin = ifelse(Type == CREs[1],
                                    -Means + StDevs, Means - StDevs),
                      ymax = ifelse(Type == CREs[1],
                                    -Means - StDevs, Means + StDevs)),
                  width = 0.5, alpha = 0.6,
                  colour = "black") +
    coord_flip() +
    scale_x_discrete(labels = kmer_labels) +
    {
      if (is.null(y_breaks)) {
        scale_y_continuous()  # Default scale without breaks
      } else {
        scale_y_continuous(breaks = y_breaks, labels = abs(y_breaks))
      }
    } +
    scale_fill_manual(values = c("turquoise", "coral"),
                      labels = CREs) +
    labs(y = "Means", x = x_label,
         title = title, fill = "CRE Type") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 15),
          text = element_text(size = rel(4.25)),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13.5),
          axis.text.y = element_text(family = "mono"),
          plot.title = element_text(size = 16, hjust = 0.5))
}

barplot_ <- function(cre_data, y_breaks = NULL, y_axis_title = "",
                     fill_legend_title = "") {
  ggplot(cre_data) +
    geom_bar(aes(x = factor(Type),
                 y = Means, fill = Type),
             stat = "identity", position = "identity",
             alpha = 0.6) +
    geom_errorbar(aes(x = factor(Type),
                      ymin = Means - StDevs,
                      ymax = Means + StDevs),
                  alpha = 0.6, width = 0.5,
                  colour = "black") +
    labs(y = y_axis_title, x = "",
         fill = fill_legend_title) +
    {
      if (is.null(y_breaks)) {
        scale_y_continuous()  # Default scale without breaks
      } else {
        scale_y_continuous(breaks = y_breaks, labels = y_breaks)
      }
    } +
    {
      if (length(unique(cre_data$Type)) > 2) {
        scale_fill_manual(values = c("Enhancer" = "turquoise",
                                     "Promoter" = "coral",
                                     "OCR" = "darkolivegreen3"))
      } else {
        scale_fill_manual(values = c("turquoise", "coral"))
      }
    } +
    theme_minimal() +
    theme(legend.position = "none",
          text = element_text(size = rel(4.5)),
          axis.title = element_text(size = 16))
}

hvioplot_ <- function(data, y_var, y_label = "",
                      fill_legend_title = "") {
  #: x_breaks <- seq(-0.12,0.12,0.06)
  ggplot(data, aes(x = 0, y = !!sym(y_var), fill = Type)) +
    geom_violinhalf(flip = 1, adjust = 0.25,
                    trim = FALSE, scale = "count",
                    position = position_dodge(width = 0),
                    linewidth = 0.1, alpha = 0.6) +
    theme_minimal() +
    #: scale_x_continuous(breaks = x_breaks ,
    #:                    labels = abs(x_breaks )) +
    scale_fill_manual(values = c("turquoise", "coral")) +
    labs(x = "", y = y_label, fill = fill_legend_title) +
    theme(legend.position = "none",
          text = element_text(size = rel(4.5)),
          axis.title = element_text(size = 13))
}

ibarplot_ <- function(cre_data, y_breaks = NULL, x_label = "",
                      y_label = "", legend_title = "") {
  ggplot(cre_data) +
    # geom_bar(aes(x = factor(Type),
    geom_bar(aes(x = Field,
                 y = Means, fill = Type),
             stat = "identity", position = position_dodge(width = 1),
             alpha = 0.6) +
    # geom_errorbar(aes(x = factor(Type),
    geom_errorbar(aes(x = Field,
                      ymin = Means - StDevs,
                      ymax = Means + StDevs, group = Type),
                  position = position_dodge(width = 1),
                  alpha = 0.6, width = 0.5,
                  colour = "black") +
    labs(y = y_label, x = x_label,
         fill = legend_title) +
    {
      if (is.null(y_breaks)) {
        scale_y_continuous()  # Default scale without breaks
      } else {
        scale_y_continuous(breaks = y_breaks, labels = y_breaks)
      }
    } +
    {
      if (length(unique(cre_data$Type)) > 2) {
        scale_fill_manual(values = c("Enhancer" = "turquoise",
                                     "Promoter" = "coral",
                                     "OCR" = "darkolivegreen3"))
      } else {
        scale_fill_manual(values = c("turquoise", "coral"))
      }
    } +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 15),
          text = element_text(size = rel(4.25)),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13.5),
          axis.text.y = element_text(family = "mono"),
          plot.title = element_text(size = 16, hjust = 0.5))
}

#______________________________________________________________________________#
#========================PRINCIPAL COMPONENTS ANALYSIS=========================#
Component <- numeric()
Variance <- numeric()
PC1 <- numeric()
PC2 <- numeric()

simple_pca <- function(data, type, n_comp = 10) {
  # Validate inputs
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data frame.")
  }
  if (length(type) != nrow(data)) {
    stop("'type' must have the same length
         as the number of rows in 'data'.")
  }

  # Scale the data
  scaled_data <- scale(data)

  # Perform PCA
  pca_result <- irlba(scaled_data, nv = n_comp)

  # Scree Plot
  sing_vals <- pca_result$d
  var_expl <- (sing_vals^2) / sum(sing_vals^2) * 100

  scree_df <- data.frame(Component = seq_along(var_expl),
                         Variance = var_expl)

  scree_plot <- ggplot(scree_df, aes(x = Component,
                                     y = Variance)) +
    geom_line(color = "cadetblue", linewidth = 0.5) +
    geom_point(color = "dodgerblue4", size = 2) +
    labs(title = "Scree Plot",
         x = "Principal Component",
         y = "% Explained Variance") +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.title = element_text(size = 13),
          axis.text = element_text(size = 11))

  # Extract scores
  scores <- pca_result$u %*% diag(pca_result$d)

  # Prepare data for ggplot
  pca_df <- data.frame(PC1 = scores[, 1],
                       PC2 = scores[, 2],
                       Type = type)

  # PCA Plot
  pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2,
                                 color = Type)) +
    geom_point(alpha = 0.6, size = 2) +
    labs(title = "PCA: First Two Components",
         x = "PC1", y = "PC2") +
    scale_color_manual(values = c("turquoise", "coral")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 15),
          text = element_text(size = rel(4.25)),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13.5),
          axis.text.y = element_text(family = "mono"),
          plot.title = element_text(size = 16, hjust = 0.5))

  # Return PCA results
  return(list(pca_result = pca_result,
              scores = scores,
              var_expl = var_expl,
              scree_plot = scree_plot,
              pca_plot = pca_plot))
}

#______________________________________________________________________________#
#===============================OUTPUT-WRAPPING================================#
count_sentwords <- function(sentences) {
  sapply(sentences, function(sentence) {
    length(strsplit(sentence, "\\s+")[[1]])
  })
}

N_digits <- function(number)
  return(str_length(as.character(number)))

wrap_output <- function(func_out, width = 50) {
  output <- capture.output(func_out)
  wrapped_output <- lapply(output, strwrap,
                           width = width)
  cat(unlist(wrapped_output), sep = "\n")
}

vectwrap <- function(func_out, width = 50, 
                     padd = 2, indexes = TRUE, round_n = FALSE) {
  form_out <- func_out
  if (is.numeric(func_out) && round_n)
    form_out <- format(round(func_out, digits = round_n), nsmall = round_n)
  wrapped <- strwrap(gsub(",", "", toString(form_out)), width = width)
  len_w <- length(wrapped)
  flag <- len_w

  if (!indexes) flag <- 0

  if (flag > 1) {
    wrapped_Ns <- count_sentwords(wrapped[1:(len_w - 1)])
    wrapped_Ns <- cumsum(c(1, wrapped_Ns))

    spc_need <- N_digits(max(wrapped_Ns)) - N_digits(wrapped_Ns)
    spc_need <- sapply(spc_need, function(x) paste(rep(" ", x), collapse = ""))

    wrapped <- paste0(spc_need, "[", wrapped_Ns, "] ", wrapped)
    wrapped[2:len_w] <- paste0("\n", wrapped[2:len_w])
  } else if (flag) {
    wrapped <- paste("[1]", wrapped)
  } else {
    wrapped[2:len_w] <- paste0("\n", wrapped[2:len_w])
  }
  #: wrapped[1] <- paste("[1]", wrapped[1])
  #: wrapped[2:len_w] <- paste("\n   ", wrapped[2:len_w])
  if (padd) {
    padding <- paste(rep("\n", padd), collapse = "")
    wrapped[1] <- paste0(padding, wrapped[1])
    wrapped[len_w] <- paste0(wrapped[len_w], padding, "\n")
  }
  cat(wrapped)
}
outwrap1 <- function(func_out, width = 50, round_n = FALSE) {
  # Capture the output of the function
  captured <- capture.output(func_out)

  # Separate names and values based on even/odd lines
  nams_capt <- captured[seq(1, length(captured), 2)] # Odd lines
  vals_capt <- captured[seq(2, length(captured), 2)] # Even lines
  vals_capt <- gsub("   ", "- -", vals_capt)
  vals_capt <- gsub("^  ", " -", vals_capt)
  vals_capt <- gsub(" $", "- ", vals_capt)

  nams_capt <- paste(nams_capt, collapse = " ")
  vals_capt <- paste(vals_capt, collapse = " ")

  # Wrap names and values independently
  wrap_nams <- strwrap(gsub(",", "", nams_capt), width = width)
  wrap_vals <- strwrap(gsub(",", "", vals_capt), width = width)
  wrap_vals <- gsub("-", " ", wrap_vals)

  # Combine wrapped names and values
  len_w <- max(length(wrap_nams), length(wrap_vals))
  wrapped <- character(len_w)
  wrapped[1] <- paste0(wrap_nams[1], "\n", wrap_vals[1])

  if (len_w > 1) {
    wrapped[2:len_w] <- paste0("\n", wrap_nams[2:len_w],
                               "\n", wrap_vals[2:len_w])
  }

  # Print the final wrapped output
  cat(wrapped, sep = "")
}

example_output <- function() {
  cat("AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC\n")
  cat(" 0   0   0   1   0   0   0   0   0   0\n")
  cat("AGG AGT ATA ATC ATG ATT\n")
  cat(" 0   0   0   0   1   0\n")
}

outwrap <- function(func_out, width = 50) {
  # Capture the output of the function
  captured <- capture.output(func_out)

  # Split the captured lines into names (odd lines) and values (even lines)
  nams_capt <- captured[seq(1, length(captured), by = 2)] # Odd lines
  vals_capt <- captured[seq(2, length(captured), by = 2)] # Even lines

  # Wrap and print the lines together
  wrapped <- character(0)
  for (i in seq_along(nams_capt)) {
    wrapped_names <- strwrap(nams_capt[i], width = width)
    wrapped_values <- strwrap(vals_capt[i], width = width)
    wrapped <- c(wrapped, wrapped_names, wrapped_values, "")
  } # Add empty line for separation

  # Print the wrapped result
  cat(wrapped, sep = "\n")
}
#______________________________________________________________________________#
#=================================TEXT-PADDING=================================#
addpadd <- function(func_out) {
  captout <- capture.output(func_out)
  captout[1] <- paste("\n", captout[1], sep = "")
  captout[2:length(captout)] <- paste("\n", captout[2:length(captout)], sep ="")
  captout <- append(captout, "\n\t")
  cat(captout)
}

hspace <- function()
    asis_output("\\textcolor{white}{\\tiny\\texttt{hi}}\\normalsize")

#______________________________________________________________________________#
#================================UNIQUE-VALUES=================================#
uniq_values <- function(vect, round_digits = 4) {
  return(round(as.numeric(names(table(vect))), digits = round_digits))
}

uniq_not <- function(vect, round_digits = 4) {
  return(round(unique(sort(vect)), digits = round_digits))
}
# ^For some reason returned duplicated values when
#  processing floating values from "gc-percentage"

#______________________________________________________________________________#
#=========================RETRIEVE-VALUES-FROM-KEYLIST=========================#
# Function to recursively collect atomic values from a nested list.
collect_atomic <- function(x) {
  if (!is.list(x)) return(x)
  result <- c()
  for (elem in x) {
    if (is.list(elem)) {
      result <- c(result, collect_atomic(elem))
    } else {
      result <- c(result, elem)
    }
  }
  return(result)
}

# Function to search for a specific key in a nested list.
search_key <- function(x, key) {
  result <- c()
  if (is.list(x)) {
    nms <- names(x)
    for (i in seq_along(x)) {
      current_name <- if (!is.null(nms)) nms[i] else ""
      # If the current element's name matches the key, collect its value(s)
      if (current_name == key) {
        if (is.list(x[[i]])) {
          result <- c(result, collect_atomic(x[[i]]))
        } else {
          result <- c(result, x[[i]])
        }
      }
      # Recursively search inside the current element (if it is a list)
      if (is.list(x[[i]])) {
        result <- c(result, search_key(x[[i]], key))
      }
    }
  }
  return(result)
}

# Function to retrieve feature values for a vector of acceptable keys.
retrieve_keylist_values <- function(nested_list, input_keys, key_order) {
  result <- c()
  for (key in input_keys) {
    result <- c(result, search_key(nested_list, key))
  }

  if (!missing(key_order))
    result <- key_order[key_order %in% result]

  return(result)
}

#______________________________________________________________________________#
#==============================RETRIEVE-SUBLIST================================#
# Function to retrieve a subset of the nested list based on acceptable keys.
retrieve_sublist <- function(x, keys) {
  # If x is not a list, there's nothing to prune.
  if (!is.list(x)) return(NULL)

  pruned <- list()
  for (nm in names(x)) {
    element <- x[[nm]]
    # If the current key is in the accepted keys, include it entirely.
    if (nm %in% keys) {
      pruned[[nm]] <- element
    } else if (is.list(element)) {
      # Otherwise, recursively prune the element.
      sub <- retrieve_sublist(element, keys)
      # Only include the element if the recursion returned a non-empty list.
      if (length(sub) > 0) {
        pruned[[nm]] <- sub
      }
    }
    # For atomic elements (non-list) not in keys, we do nothing.
  }
  pruned
}

# Function to recursively retrieve all the leaves (atomic values)
# from a nested list.
retrieve_leaves <- function(x) {
  # If x is a list, recurse into its elements.
  if (is.list(x)) {
    leaves <- c()
    for (item in x) {
      leaves <- c(leaves, retrieve_leaves(item))
    }
    leaves
  } else {
    # x is atomic (or an atomic vector) so return it.
    x
  }
}

retrieve_leaves_with_names <- function(x, current_name = NULL) {
  if (!is.list(x)) {
    # x is atomic; if a current name exists, assign it to all elements
    n <- length(x)
    if (!is.null(current_name)) {
      names(x) <- rep(current_name, n)
    }
    return(x)
  }

  result <- c()
  nms <- names(x)

  for (i in seq_along(x)) {
    # If the current element has a name, use it; otherwise, use the parent's name.
    key <- if (!is.null(nms) && nms[i] != "") nms[i] else current_name
    result <- c(result, retrieve_leaves_with_names(x[[i]], key))
  }

  result
}
#______________________________________________________________________________#
#==============================PRINT-LIST-AS-TREE==============================#
print_tree <- function(x, prefix = "") {
  if (is.list(x)) {
    n <- length(x)
    keys <- names(x)
    for (i in seq_along(x)) {
      is_last <- i == n
      # Choose the branch symbol depending on position.
      branch <- if (is_last) "└── " else "├── "

      # Check if the current element has a name.
      current_name <- if (!is.null(keys) && keys[i] != "") keys[i] else NULL

      if (!is.null(current_name)) {
        # Print the name of the current branch.
        cat(prefix, branch, current_name, "\n", sep = "")
        # Update the prefix for children.
        new_prefix <- paste0(prefix, if (is_last) "    " else "│   ")
        child <- x[[i]]
        # If the child is a list, recurse.
        if (is.list(child)) {
          print_tree(child, prefix = new_prefix)
        } else if (is.atomic(child) && length(child) > 1) {
          # If the child is an atomic vector with multiple elements,
          # print each element on a new line.
          m <- length(child)
          for (j in seq_along(child)) {
            leaf_is_last <- j == m
            leaf_branch <- if (leaf_is_last) "└── " else "├── "
            cat(new_prefix, leaf_branch, child[j], "\n", sep = "")
          }
        } else {
          # For a single atomic value, print it as a leaf.
          cat(new_prefix, "└── ", child, "\n", sep = "")
        }
      } else {
        # For unnamed elements (like "MT" or "SEN"), print directly.
        cat(prefix, branch, x[[i]], "\n", sep = "")
      }
    }
  } else if (is.atomic(x) && length(x) > 1) {
    # If x is an atomic vector (e.g. a character vector with many values),
    # print each element.
    m <- length(x)
    for (i in seq_along(x)) {
      is_last <- i == m
      branch <- if (is_last) "└── " else "├── "
      cat(prefix, branch, x[i], "\n", sep = "")
    }
  } else {
    # Otherwise, simply print the value.
    cat(prefix, x, "\n", sep = "")
  }
}

print_tree_names <- function(x, prefix = "",
                             file_name = NULL, file_conn = NULL) {
  # Check if a file name is provided to create a new file
  if (!is.null(file_name)) {
    if (!is.null(file_conn)) {
      stop("Cannot provide both 'file_name' and 'file_conn'")
    }
    # Open a connection to the file in write mode (creates or overwrites)
    file_conn <- file(file_name, "w")
    # Ensure the connection is closed after the function exits
    on.exit(close(file_conn), add = TRUE)
  }

  if (is.list(x)) {
    n <- length(x)
    nms <- names(x)
    for (i in seq_along(x)) {
      is_last <- i == n
      # A leaf is defined as an element that is not a list.
      is_leaf <- !is.list(x[[i]])

      # Use different branch symbols for leaves.
      branch <- if (is_leaf) {
        if (is_last) "└──=> " else "├──=> "
      } else {
        if (is_last) "└── " else "├── "
      }

      # Retrieve the name if available.
      current_name <- if (!is.null(nms) && nms[i] != "") nms[i] else NULL

      if (!is.null(current_name)) {
        if (!is.null(file_conn)) {
          cat(prefix, branch, current_name, "\n", sep = "", file = file_conn)
        } else {
          cat(prefix, branch, current_name, "\n", sep = "")
        }
      }

      # Only recurse if the element is a list.
      if (is.list(x[[i]])) {
        new_prefix <- paste0(prefix, if (is_last) "    " else "│   ")
        print_tree_names(x[[i]], new_prefix, file_conn = file_conn)
      }
    }
  }
}

example_tree <-
  list(WS = list(basepercs = "A",
                 temp = "B",
                 shan = "C",
                 SA = list(RA = list(r_id = "D",
                                     r_sc = "E"),
                           RCA = list(rc_id = "F",
                                      rc_sc = "G"))),
       KS = list(KD = list(kx_shan = "H",
                           kx_adiv = "I",
                           kx_rdiv = "J")),
       EK = list(ki_prod = "K",
                 ki_bcds = "L",
                 ki_bclp = "M"),
       RK = list(dna_shape = "N",
                 RKS = list(ki_sp_pol = "O",
                            ki_sp_lag = "P")))

#______________________________________________________________________________#
#===========================SUBSTITUTE-NAMES-IN-TREE===========================#
substitute_names <- function(x, mapping) {
  # If x is not a list, nothing to do.
  if (!is.list(x)) return(x)

  # If the list has names, substitute them using the mapping.
  if (!is.null(names(x))) {
    names(x) <- sapply(names(x), function(nm) {
      if(nm %in% names(mapping)) mapping[[nm]] else nm
    })
  }

  # Apply recursively to each element.
  for (i in seq_along(x)) {
    if (is.list(x[[i]])) {
      x[[i]] <- substitute_names(x[[i]], mapping)
    }
  }
  return(x)
}
