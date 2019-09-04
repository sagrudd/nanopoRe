

#' prepare an 3x infographic-styled summary plot
#'
#' This method will prepare a 3x graphical plot summarising user data with an emphasis on visualisation
#'
#' @importFrom emojifont fontawesome
#' @import emojifont
#' @usage infoGraphicPlot3(identifier, panelA, panelB, panelC,
#'     reportDPI = 90)
#' @param identifier a key for the plot
#' @param panelA information for panel1 in format of panelA=c(key='', value='', icon='')
#' @param panelB see panelA
#' @param panelC see panelA
#' @param reportDPI resolution for plot
#' @return link to a figure for display
#'
#' @examples
#' library(emojifont)
#' infoGraphicPlot3(identifier='mappingSummary',
#'     panelA=c(key='Alignments', value='123456', icon='fa-sliders'),
#'     panelB=c(key='Mapping Yield', value='680M', icon='fa-calculator'),
#'     panelC=c(key='Average Accuracy', value='89.0 %', icon='fa-area-chart'))
#'
#' @export
infoGraphicPlot3 <- function(identifier, panelA, panelB, panelC, reportDPI = 90) {
    # panelA=c(key='', value='', icon='')
    figures <- 3

    df <- data.frame(x = cumsum(c(2, rep(6.5, figures - 1))), y = rep(2, figures), h = rep(4, figures),
        w = rep(6, figures))

    df$key <- c(panelA[["key"]], panelB[["key"]], panelC[["key"]])
    df$info <- c(panelA[["value"]], panelB[["value"]], panelC[["value"]])
    df$icon <- emojifont::fontawesome(c(panelA[["icon"]], panelB[["icon"]], panelC[["icon"]]))

    df$colour <- rep("steelblue", figures)

    InfoGraphicValueBoxes <- ggplot(df, aes_string("x", "y", height = "h", width = "w", label = "key",
        fill = "colour")) + geom_tile(fill = brewer.pal(9, "Blues")[7]) + geom_text(color = brewer.pal(9,
        "Blues")[3], hjust = "left", nudge_y = -1.5, nudge_x = -2.6, size = 5) + geom_text(label = df$icon,
        family = "fontawesome-webfont", colour = brewer.pal(9, "Blues")[5], size = 23, hjust = "right",
        nudge_x = 2.85, nudge_y = 0.8) + geom_text(label = df$info, size = 10, color = brewer.pal(9,
        "Blues")[2], fontface = "bold", nudge_x = -2.6, hjust = "left") + coord_fixed() + scale_fill_brewer(type = "qual",
        palette = "Dark2") + theme_void() + guides(fill = FALSE)

    infographicFile <- file.path(getRpath(), paste0(identifier, ".igraphic.png"))

    ggplot2::ggsave(infographicFile, plot = InfoGraphicValueBoxes, device = "png", units = "cm", width = 25,
        height = 5, dpi = reportDPI)

    # knitr::include_graphics(infographicFile)
    return(infographicFile)
}




#' prepare a 4x infographic-styled summary plot
#'
#' This method will prepare a 4x graphical plot summarising user data with an emphasis on visualisation
#'
#' @importFrom emojifont fontawesome
#' @import emojifont
#' @usage infoGraphicPlot4(identifier, panelA, panelB, panelC,
#'     panelD, reportDPI = 90)
#' @param identifier a key for the plot
#' @param panelA information for panel1 in format of panelA=c(key='', value='', icon='')
#' @param panelB see panelA
#' @param panelC see panelA
#' @param panelD see panelA
#' @param reportDPI resolution for plot
#' @return link to a figure for display
#'
#' @examples
#' library(emojifont)
#' infoGraphicPlot4(identifier='mappingSummary',
#'     panelA=c(key='Alignments', value='123456', icon='fa-sliders'),
#'     panelB=c(key='Mapping Yield', value='680M', icon='fa-calculator'),
#'     panelC=c(key='Average Accuracy', value='89.0 %', icon='fa-area-chart'),
#'     panelD=c(key='Average Identity', value='96.2 %', icon='fa-info'))
#'
#' @export
infoGraphicPlot4 <- function(identifier, panelA, panelB, panelC, panelD, reportDPI = 90) {
    # panelA=c(key='', value='', icon='')
    figures <- 4

    df <- data.frame(x = cumsum(c(2, rep(6.5, figures - 1))), y = rep(2, figures), h = rep(4, figures),
        w = rep(6, figures))

    df$key <- c(panelA[["key"]], panelB[["key"]], panelC[["key"]], panelD[["key"]])
    df$info <- c(panelA[["value"]], panelB[["value"]], panelC[["value"]], panelD[["value"]])
    df$icon <- emojifont::fontawesome(c(panelA[["icon"]], panelB[["icon"]], panelC[["icon"]], panelD[["icon"]]))

    df$colour <- rep("steelblue", figures)

    InfoGraphicValueBoxes <- ggplot(df, aes_string("x", "y", height = "h", width = "w", label = "key",
        fill = "colour")) + geom_tile(fill = brewer.pal(9, "Blues")[7]) + geom_text(color = brewer.pal(9,
        "Blues")[3], hjust = "left", nudge_y = -1.5, nudge_x = -2.6, size = 3.5) + geom_text(label = df$icon,
        family = "fontawesome-webfont", colour = brewer.pal(9, "Blues")[5], size = 14, hjust = "right",
        nudge_x = 2.85, nudge_y = 0.9) + geom_text(label = df$info, size = 9, color = brewer.pal(9,
        "Blues")[2], fontface = "bold", nudge_x = -2.6, hjust = "left") + coord_fixed() + scale_fill_brewer(type = "qual",
        palette = "Dark2") + theme_void() + guides(fill = FALSE)

    infographicFile <- file.path(getRpath(), paste0(identifier, ".igraphic.png"))

    ggplot2::ggsave(infographicFile, plot = InfoGraphicValueBoxes, device = "png", units = "cm", width = 25,
        height = 5, dpi = reportDPI)

    # knitr::include_graphics(infographicFile)
    return(infographicFile)
}




#' prepare a 5x infographic-styled summary plot
#'
#' This method will prepare a 5x graphical plot summarising user data with an emphasis on visualisation
#'
#' @importFrom emojifont fontawesome
#' @import emojifont
#' @usage infoGraphicPlot5(identifier, panelA, panelB, panelC,
#'     panelD, panelE, reportDPI = 90)
#' @param identifier a key for the plot
#' @param panelA information for panel1 in format of panelA=c(key='', value='', icon='')
#' @param panelB see panelA
#' @param panelC see panelA
#' @param panelD see panelA
#' @param panelE see panelA
#' @param reportDPI resolution for plot
#' @return link to a figure for display
#'
#' @examples
#' library(emojifont)
#' infoGraphicPlot4(identifier='mappingSummary',
#'     panelA=c(key='Alignments', value='123456', icon='fa-sliders'),
#'     panelB=c(key='Mapping Yield', value='680M', icon='fa-calculator'),
#'     panelC=c(key='Average Accuracy', value='89.0 %', icon='fa-area-chart'),
#'     panelD=c(key='Average Identity', value='96.2 %', icon='fa-info'),
#'     panelE=c(key='Average Identity', value='96.2 %', icon='fa-info'))
#'
#' @export
infoGraphicPlot5 <- function(identifier, panelA, panelB, panelC, panelD, panelE, reportDPI = 90) {
    # panelA=c(key='', value='', icon='')
    figures <- 5

    df <- data.frame(x = cumsum(c(2, rep(6.5, figures - 1))), y = rep(2, figures), h = rep(4, figures),
        w = rep(6, figures))

    df$key <- c(panelA[["key"]], panelB[["key"]], panelC[["key"]], panelD[["key"]], panelE[["key"]])
    df$info <- c(panelA[["value"]], panelB[["value"]], panelC[["value"]], panelD[["value"]], panelE[["value"]])
    df$icon <- emojifont::fontawesome(c(panelA[["icon"]], panelB[["icon"]], panelC[["icon"]], panelD[["icon"]],
        panelE[["icon"]]))

    df$colour <- rep("steelblue", figures)

    InfoGraphicValueBoxes <- ggplot(df, aes_string("x", "y", height = "h", width = "w", label = "key",
        fill = "colour")) + geom_tile(fill = brewer.pal(9, "Blues")[7]) + geom_text(color = brewer.pal(9,
        "Blues")[3], hjust = "left", nudge_y = -1.5, nudge_x = -2.6, size = 3.5) + geom_text(label = df$icon,
        family = "fontawesome-webfont", colour = brewer.pal(9, "Blues")[5], size = 13.3, hjust = "right",
        nudge_x = 2.85, nudge_y = 0.8) + geom_text(label = df$info, size = 9, color = brewer.pal(9,
        "Blues")[2], fontface = "bold", nudge_x = -2.6, hjust = "left") + coord_fixed() + scale_fill_brewer(type = "qual",
        palette = "Dark2") + theme_void() + guides(fill = FALSE)

    infographicFile <- file.path(getRpath(), paste0(identifier, ".igraphic.png"))

    ggplot2::ggsave(infographicFile, plot = InfoGraphicValueBoxes, device = "png", units = "cm", width = 25,
        height = 5, dpi = reportDPI)

    # knitr::include_graphics(infographicFile)
    return(infographicFile)
}

