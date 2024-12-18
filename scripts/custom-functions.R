library(knitr)
library(kableExtra)
library(stringi)
library(dplyr)

library(ggplot2)
library(cowplot)
library(see)

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

# To plot pyramid-plots
pyrplot_ <- function(cre_data, kmer_labels, title,
                     x_label, y_label, y_breaks) {
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
    scale_y_continuous(breaks = y_breaks,
                       labels = abs(y_breaks)) +
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

barplot_ <- function(cre_data, y_breaks, y_axis_title = "",
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
    scale_y_continuous(breaks = y_breaks,
                       labels = y_breaks) +
    scale_fill_manual(values = c("turquoise", "coral")) +
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

#______________________________________________________________________________#
#===============================OUTPUT-WRAPPING================================#
wrap_output <- function(func_out, width = 50) {
  output <- capture.output(func_out)
  wrapped_output <- lapply(output, strwrap,
                           width = width)
  cat(wrapped_output, sep = "\n")
}

outputwrap <- function(func_out, width = 50) {
  form_out <- format(round(func_out, digits = 2),
                     nsmall = 2)
  wrapped <- strwrap(gsub(",", "", toString(form_out)),
                     width = width)
  wrapped[1] <- paste("\n\n[1]", wrapped[1])
  len_w <- length(wrapped)
  wrapped[2:len_w] <- paste("\n   ", wrapped[2:len_w])
  cat(wrapped)
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



