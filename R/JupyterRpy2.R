


#' method to mask R output by writing to a log file
#'
#' creates a log file and writes R code and messages on STDERR to the log file
#'
#' @param filename the path to the logfile to use
#' @param append is a boolean defining whether the log file should be appended
#' or ab initio
#' @return None
#'
#' @examples
#' setLogFile('demofile.log', append=FALSE)
#' unsetLog()
#'
#' @export
setLogFile <- function(filename = NULL, append = TRUE) {
    if (is.null(filename)) {
        if (hasCachedObject("logfile")) {
            filename <- getCachedObject("logfile")
        } else {
            filename = "unnamed.log"
        }
    } else {
        setCachedObject("logfile", filename)
    }
    con <- file(filename)
    sink(con, append = append)
    sink(con, append = append, type = "message")
}

#' method to output capture to log file
#'
#' this assumes that a log file has been written to and ends the sinking
#'
#' @return None
#'
#' @examples
#' setLogFile('demofile.log', append=FALSE)
#' unsetLog()
#'
#' @export
unsetLog <- function() {
    sink()
    sink(type = "message")
}
