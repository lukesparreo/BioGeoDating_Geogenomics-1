# Example data
#newick_string <- "((A[&index=1,a=1]:1,B[&index=2,a=2]:1)[&index=3,a=3]:1,C[&index=4,a=4]:1);"
newick_string <- ase@treetext

# Example tibble
library(tibble)
#tibble_data <- tibble(index = 1:4, height_0.95_HPD = list(c(0.1, 0.1), c(0.2, 0.2), c(0.3, 0.3), c(0.4, 0.4)))

# Real data 
tibble_data <- time_tree_data[c(1,3)]

# Load the stringr package
library(stringr)

# Function to replace index values with corresponding height values
replace_height <- function(match) {
  index <- as.numeric(gsub("[^0-9]", "", match))  # Extract index value from the match
  
  # Check if index exists in tibble
  if (index %in% tibble_data$index) {
    hpd_values <- tibble_data$age_0.95_HPD[tibble_data$index == index][[1]]
    paste0("&index=", index, ",height_0.95_HPD={", paste(hpd_values, collapse = ","), "}")
  } else {
    warning(paste("Index", index, "not found in tibble. Replacement skipped."))
    match  # If index not found, return the original match
  }
}

# Replace index values in the newick string
modified_string <- str_replace_all(newick_string, "&index=\\d+", replace_height)

# Print the modified string
print(modified_string)


