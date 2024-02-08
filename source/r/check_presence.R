library(dplyr)
library(purrr)

# Set working directory
setwd("YOUR_WORKING_DIRECTORY")

# West-Europa
# België, Frankrijk, Duitsland, Luxemburg, Nederland, Zwitserland, Oostenrijk


check_presence_from_file <- function(input_file = "./EJP_inse01_species_unique.txt", countries = c("BE", "FR", "DE", "LU", "NL", "CH", "AT")) {
  # Read scientific names from the input file
  scientificNames <- readLines(input_file)

  # Remove any empty lines
  scientificNames <- scientificNames[scientificNames != ""]

  # Initialize vector to store presence results
  presence <- logical(length(scientificNames))

  # Loop over each country
  for (country in countries) {
    # Function to check occurrence data and handle errors
    check_occurrence <- function(name) {
      data <- tryCatch({
        rgbif::occ_data(
          scientificName = name,
          country = country,  # Specify one country at a time
          limit = 1
        )
      }, error = function(e) {
        return(list(error = TRUE))
      })

      if (inherits(data, "error") || is.null(data$data)) {
        return(FALSE)  # Occurrence data is not present or error occurred
      } else {
        return(TRUE)   # Occurrence data is present
      }
    }

    # Map over scientific names and check occurrence data for the current country
    country_presence <- purrr::map_lgl(scientificNames, check_occurrence)

    # Update presence vector with results for the current country
    presence <- presence | country_presence
  }

  # Create tibble with results
  result_table <- tibble::tibble(present = presence, scientificName = scientificNames)

  return(result_table)
}

# Example usage:
result_WestEurope <- check_presence_from_file()
print(result_WestEurope)

write.csv(result_WestEurope, "result_WestEurope.csv")
