% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InfographicPlot.R
\name{infoGraphicPlot4}
\alias{infoGraphicPlot4}
\title{prepare a 4x infographic-styled summary plot}
\usage{
infoGraphicPlot4(identifier, panelA, panelB, panelC,
    panelD, reportDPI = 90)
}
\arguments{
\item{identifier}{a key for the plot}

\item{panelA}{information for panel1 in format of panelA=c(key='', value='',
icon='')}

\item{panelB}{see panelA}

\item{panelC}{see panelA}

\item{panelD}{see panelA}

\item{reportDPI}{resolution for plot}
}
\value{
link to a figure for display
}
\description{
This method will prepare a 4x graphical plot summarising user data with an
emphasis on visualisation
}
\examples{
library(emojifont)
infoGraphicPlot4(identifier='mappingSummary',
    panelA=c(key='Alignments', value='123456', icon='fa-sliders'),
    panelB=c(key='Mapping Yield', value='680M', icon='fa-calculator'),
    panelC=c(key='Average Accuracy', value='89.0 \%', icon='fa-area-chart'),
    panelD=c(key='Average Identity', value='96.2 \%', icon='fa-info'))

}
