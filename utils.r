#' get_requirements installs the packages from the given file.
#'
#' @param file the requirements file.
get_requirements <- function(file) {
  if (Sys.getenv("PACKAGES_INSTALLED") != "") {
    return()
  }

  if (!require("requiRements")) {
    install.packages("requiRements")
  }
  requiRements::install(NULL, file)

  Sys.setenv(PACKAGES_INSTALLED = TRUE)
}
