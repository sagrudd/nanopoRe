
#' set path for R results
#'
#' This method sets a path where default R results will be stored during an
#' analysis
#'
#' @param path is the location to the directory for storing results
#' @return None
#'
#' @examples
#' setRpath(file.path(tempdir(), "nanopoRe"))
#' getRpath()
#'
#' @export
setRpath <- function(path) {
    assign("r_path", path, envir = get(getEnvironment()))
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
}



#' get path for R results
#'
#' This method returns the path where default R results will be stored
#' during an analysis
#'
#' @return character representation of path
#'
#' @examples
#' setRpath(file.path(tempdir(), "nanopoRe"))
#' getRpath()
#'
#' @export
getRpath <- function() {
    return(get("r_path", envir = get(getEnvironment())))
}



#' get the defined output format for ggplot2 figures
#'
#' This method defines the expected output format for figure management
#'
#' @return character of format
#'
#' @examples
#' getOutputFormat()
#' setOutputFormat(listOutputFormat()[2])
#' getOutputFormat()
#'
#' @export
getOutputFormat <- function() {
    objName = "gg2out"
    if (hasCachedObject(objName)) {
        return(getCachedObject(objName))
    } else {
        return(listOutputFormat()[1])
    }
}

#' set the defined output format for ggplot2 figures
#'
#' This method defines the expected output format for figure management
#'
#' @param gg2out string specifying type
#' @return character of format
#'
#' @examples
#' getOutputFormat()
#' setOutputFormat(listOutputFormat()[2])
#' getOutputFormat()
#'
#' @export
setOutputFormat <- function(gg2out) {
    objName = "gg2out"
    if (gg2out %in% listOutputFormat()) {
        setCachedObject(objName, gg2out)
    } else {
        warning(paste0(gg2out,
            "is not a valid outputformat - see listOutputFormat()"))
    }
}


#' list available output frameworks for ggplot2 figures
#'
#' This method presents the list of allowed data types for managing the
#' simplified presentation of ggplot2 figures
#'
#' @return vector of allowed types
#'
#' @examples
#' listOutputFormat()
#'
#' @export
listOutputFormat <- function() {
    return(c("raw", "file", "knitr", "jupyter"))
}



ggplot2save <- function(plot, filename=NA) {
    dim <- getPlotDimensions()
    dest <- tempfile(pattern = "", tmpdir = getRpath(), fileext = ".png")
    if (!is.na(filename)) {
        dest <- file.path(getRpath(), filename)
    }
    ggplot2::ggsave(
        dest, plot = plot, width = dim$width, height = dim$height,
        units = dim$units, dpi = dim$dpi)
    return(dest)
}

#' @importFrom IRdisplay display_png
ggplot2handler <- function(plot) {
    if (getOutputFormat() == "raw") {
        return(plot)
    } else if (getOutputFormat() == "file") {
        return(ggplot2save(plot))
    } else if (getOutputFormat() == "knitr") {
        return(plot)
    } else if (getOutputFormat() == "jupyter") {
        return(IRdisplay::display_png(file = ggplot2save(plot)))
    }
    message(paste0(
        "ggplot2 output format { ",
        getOutputFormat(), "} not known - returning *raw*\n"))
    return(plot)
}


#' return defined dimensions for an allowed ggplot2 object
#'
#' This method presents the list with the defined metrics for a
#' ggplot2 presentation
#'
#' @return list defining dpi, width, height and units
#'
#' @examples
#' getPlotDimensions()
#'
#' @export
getPlotDimensions <- function() {
    objName = "gg2dimensions"
    if (hasCachedObject(objName)) {
        return(getCachedObject(objName))
    } else {
        return(list(dpi = 90, width = 9, height = 6, units = "in"))
    }
}


#' set dimension definitions for an allowed ggplot2 object
#'
#' This method sets and stores the list with the defined metrics for a
#' ggplot2 presentation
#'
#' @param dim list of properties
#' @return list defining dpi, width, height and units
#'
#' @examples
#' setPlotDimensions(list(dpi=180, width=21, height=12, units='cm'))
#' getPlotDimensions()
#'
#' @export
setPlotDimensions <- function(dim) {
    if (is.list(dim)) {
        objName = "gg2dimensions"
        keys <- c("dpi", "width", "height", "units")
        if (length(which(!names(dim) %in% keys)) > 0) {
            dim[which(!names(dim) %in% keys)] <- NULL
        }
        if (length(which(!keys %in% names(dim))) > 0) {
            warning(paste0(
                "Cannot update plot dimensions [",
                keys[which(!keys %in% names(dim))], "] missing"))
            return()
        }

        if (is.numeric(dim$dpi) && is.numeric(dim$width) &&
            is.numeric(dim$height) && is.character(dim$units)) {
            setCachedObject(objName, dim)
        } else {
            warning(paste0("Cannot update plot dimensions - ",
                "one or more keys of wrong type"))
        }
    }
}
