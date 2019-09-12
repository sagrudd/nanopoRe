

#' get git metadata for current project
#'
#' The ONT tutorials are typically checked out of github - there is a potential
#' for misalignment of a user's working branch and current - would be sensible
#' to include a brief git comment so that Techical Support can follow-up with
#' assessment of versioning. This code will not work with code that has been
#' extracted from either CONDA or Bioconductor ...
#'
#' @return kable
#'
#' @importFrom git2r in_repository
#' @importFrom git2r branches
#' @importFrom git2r status
#' @examples
#' init()
#' getGitDescription()
#'
#' @export
getGitDescription <- function() {
    if (git2r::in_repository()) {
        git2r::branches()
        git2r::status()
    }
}
