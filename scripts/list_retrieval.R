# =========================================================================== #
#   RETRIEVE KEYLIST VALUES: From key set (names), retrieves values			  #
#                            associated (inside a named list of lists).		  #
# =========================================================================== #

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
retireve_values <- function(nested_list, input_keys, key_order) {
  result <- c()
  for (key in input_keys) {
    result <- c(result, search_key(nested_list, key))
  }

  if (!missing(key_order))
    result <- key_order[key_order %in% result]

  return(result)
}

# Example nested list.
feat_func <- list(
  WS = list(
    basepercs = "A",
    temp = "B",
    shan = "C",
    SA = list(
      RA = list(r_id = "D", r_sc = "E"),
      RCA = list(rc_id = "F", rc_sc = "G")
    )
  ),
  KS = list(
    KD = list(kx_shan = "H", kx_adiv = "I", kx_rdiv = "J")
  ),
  EK = list(
    ki_prod = "K",
    ki_bcds = "L",
    ki_bclp = "M"
  )
)

# Curated input vector.
input_keys <- c("basepercs", "rc_id", "KS", "ki_bcds")

# Retrieve the corresponding values.
result <- retrieve_feature_values(feat_func, input_keys)
print(result)


# =========================================================================== #
#   RETRIEVE TREE FROM LEAVES: Assuming a "tree" is a list of lists, this	  #
#							   function	retrieves only the branches			  #
#							   associated to a given set of "leaves" (values).#
# =========================================================================== #

# Function to retrieve a subset of the nested list that only contains branches
# leading to allowed leaves.
retrieve_sublist <- function(x, allowed) {
  # If x is atomic (or not a list), filter to allowed elements.
  if (!is.list(x)) {
    allowed_elems <- x[x %in% allowed]
    if (length(allowed_elems) > 0) return(allowed_elems) else return(NULL)
  }

  new_list <- list()

  # Iterate over elements of the list (using names if available)
  if (is.null(names(x))) {
    for (i in seq_along(x)) {
      res <- retrieve_sublist(x[[i]], allowed)
      if (!is.null(res)) {
        new_list[[i]] <- res
      }
    }
  } else {
    for (nm in names(x)) {
      res <- retrieve_sublist(x[[nm]], allowed)
      if (!is.null(res)) {
        new_list[[nm]] <- res
      }
    }
  }

  if (length(new_list) > 0) return(new_list) else return(NULL)
}

# Function to retrieve all the leaves (atomic values) from a nested list.
retrieve_leaves <- function(x) {
  if (!is.list(x)) return(x)
  leaves <- c()
  for (item in x) {
    leaves <- c(leaves, retrieve_leaves(item))
  }
  leaves
}

# =========================================================================== #
#   PRINT LIST AS TREE: Assuming a "tree" is a list of lists, this function   #
#						function prints the structure in a hierarchycal way.  #
# =========================================================================== #

print_listtree <- function(x, prefix = "") {
  if (is.list(x)) {
    n <- length(x)
    nms <- names(x)
    for (i in seq_along(x)) {
      is_last <- i == n
      branch <- if (is_last) "└── " else "├── "

      # Retrieve the name if available; otherwise, set to NULL.
      current_name <- if (!is.null(nms) && nms[i] != "") nms[i] else NULL

      # Print the name only if it exists.
      if (!is.null(current_name)) {
        cat(prefix, branch, current_name, "\n", sep = "")
      }

      # Prepare new prefix for children.
      new_prefix <- paste0(prefix, if (is_last) "    " else "│   ")

      # Recurse if the element is a list.
      if (is.list(x[[i]])) {
        print_listtree(x[[i]], new_prefix)
      }
    }
  }
}
