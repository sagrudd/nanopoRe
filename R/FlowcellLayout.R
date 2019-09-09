

SequencingSummaryGetChannelMap <- function(platform) {
    if (platform == "MinION") {
        return(getMinIONChannelMap())
    } else if (platform == "Flongle") {
        return(getFlongleChannelMap())
    } else if (platform == "PromethION") {
        return(getPromethIONChannelMap())
    }
    return(NULL)
}



#' produce the channelMap for a MinION flowcell for spatial plots
#'
#' prepares a data.frame of channelIds and their X, Y coordinates for a MinION flowcell
#'
#' @return data.frame with channel, row and col columns
#'
#' @examples
#' head(getMinIONChannelMap())
#'
#' @export
getMinIONChannelMap <- function() {
    # build the map for R9.4.1 flowcell, as a long-form dataframe
    blockCalc <- function(i) {
        m <- matrix(seq(i, i + 63, by = 1), ncol = 8, byrow = TRUE)
        cbind(m[seq(5, 8, by = 1), ], m[seq(4), rev(seq(8))])
    }
    layout <- do.call(rbind, lapply(c(1, 449, 385, 321, 257, 193, 129, 65), blockCalc))
    channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)), which(layout == as.vector(layout),
        arr.ind = TRUE)))
    return(channelMap)
}


#' produce the channelMap for a flongle flowcell for spatial plots
#'
#' prepares a data.frame of channelIds and their X, Y coordinates for a flongle flowcell
#'
#' @return data.frame with channel, row and col columns
#'
#' @examples
#' head(getFlongleChannelMap())
#'
#' @export
getFlongleChannelMap <- function() {
    layout <- matrix(c(seq(1, 12), 0, seq(13, 24), 0, seq(25, 114), 0, seq(115, 126), 0), ncol = 13,
        byrow = TRUE)
    layout <- layout[rev(seq(10)), ]
    channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)), which(layout == as.vector(layout),
        arr.ind = TRUE)))
    return(channelMap)
}


#' produce the channelMap for a PromethION flowcell for spatial plots
#'
#' prepares a data.frame of channelIds and their X, Y coordinates for a PromethION flowcell
#'
#' @return data.frame with channel, row and col columns
#'
#' @examples
#' head(getPromethIONChannelMap())
#'
#' @export
getPromethIONChannelMap <- function() {
    chunk <- function(i) {
        m <- matrix(seq_len(250), ncol=10, byrow=TRUE)
        m + i
    }
    layout <- do.call(cbind, lapply(seq(from=0, to=2750, by=250), chunk))
    channelMap <- as.data.frame(cbind(channel = as.vector(t(layout)), which(layout == as.vector(layout),
                                                                            arr.ind = TRUE)))
    return(channelMap)
}

