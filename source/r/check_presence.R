#' Deze functie gaat na of een soort volgens GBIF voorkomt in een gekozen set
#' van landen # nolint start
#'
#' Default zijn Europese landen ingesteld: Albania, Andorra, Austria, Belarus, Belgium, Bosnia and Herzegovina, Bulgaria, Croatia, Czechia, Denmark, Estonia, Finland, France, Germany, Greece, Hungary, Ireland, Italy, Kosovo, Latvia, Liechtenstein, Lithuania, Luxembourg, Monaco, Montenegro, Netherlands, North Macedonia, Norway, Poland, Portugal, Romania, San Marino, Serbia, Slovakia, Slovenia, Spain, Sweden, Switzerland, Ukraine, United Kingdom
#'
#' @param input A character vector of Scientific names or a .txt file from which
#' Scientific Names can be read.
#' @param countries a character vector of country codes
#'
#' @example
#' result_Europe <- check_presence_from_file()
#' print(result_Europe) # nolint end
check_presence <- function(
    input = NULL,
    countries = c("AL", "AD", "AT", "BY", "BE", "BA", "BG", "HR", "CZ", "DK", "EE", "FI", "FR", "DE", "GR"  ,"HU", "IE", "IT", "XK", "LV", "LI", "LT", "LU", "MC", "ME", "NL", "MK", "NO", "PL", "PT", "RO", "SM", "RS", "SK", "SI", "ES", "SE", "CH", "UA", "UK")) {
  require(dplyr)
  require(purrr)
  if (length(input) == 1 && !is.null(input) && file.exists(input)) {
    # Read scientific names from the input file
    scientific_names <- readLines(input)
  } else if (is.character(input)) {
    scientific_names <- input
  } else {
    stop("Input must be either a character vector or a valid file path.")
  }

  # Remove any empty lines
  scientific_names <- scientific_names[scientific_names != ""]

  # Initialize vector to store presence results
  presence <- logical(length(scientific_names))

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

    # Map over scientific names and check occurrence data for the current
    # country
    country_presence <- purrr::map_lgl(scientific_names, check_occurrence)

    # Update presence vector with results for the current country
    presence <- presence | country_presence
  }

  # Create tibble with results
  result_table <- tibble::tibble(
    present = presence,
    scientific_name = scientific_names
  )

  return(result_table)
}
