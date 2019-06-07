
#' set path for reference genome object
#'
#' This method defines the path where reference genome is stored
#'
#' @param reference.file is the location of the fasta format sequence file for reference genome
#' @return None
#'
#' @examples
#' \dontrun{
#' setReferenceGenome(file.path("ReferenceData", "human_g1k_v37.fasta"))
#' }
#'
#' @export
setReferenceGenome <- function(reference.file) {
  assign("reference.file", reference.file, envir=get(getEnvironment()))
}

#' get path reference genome object
#'
#' This method returns the stored path for the reference genome object
#'
#' @return character representation of path to stored fasta format file object
#'
#' @examples
#' \dontrun{
#' getReferenceGenome()
#' }
#'
#' @export
getReferenceGenome <- function() {
  return(get("reference.file", envir=get(getEnvironment())))
}



#' load reference genome into memory
#'
#' This method will parse the defined reference.file into a DNAStringSet object for handling in memory
#'
#' @return NULL
#'
#' @examples
#' \dontrun{
#' getReferenceGenome()
#' }
#'
#' @export
loadReferenceGenome <- function() {
  message(paste0("loading reference genome(",getReferenceGenome(),")"))
  # derive a referenceGenome object from the named fasta elements in the provided fasta reference resource
  referenceGenomeSequence <- readDNAStringSet(getReferenceGenome())
  referenceGenome <- data.frame(id=gsub(" .+", "", names(referenceGenomeSequence)),sid=seq_along(names(referenceGenomeSequence)),stringsAsFactors = FALSE)
  referenceGenome$sid <- seq(nrow(referenceGenome))
  assign("referenceGenome", referenceGenome, envir=get(getEnvironment()))
  assign("referenceGenomeSequence", referenceGenomeSequence, envir=get(getEnvironment()))
}


#' cleanup the reference genome data
#'
#' This method will remove unwanted chromosomes from the reference genome collection; useful for when analyses are to
#' be linked to a set of core autosomes, for example.
#'
#' @param delIds is a vector of chromosomes to be dropped
#' @return None
#'
#' @examples
#' \dontrun{
#' cleanReferenceGenome(chrId)
#' }
#'
#' @export
cleanReferenceGenome <- function(delIds) {
  referenceGenome <- get("referenceGenome", envir=get(getEnvironment()))
  overlap <- match(delIds, referenceGenome[,1])
  if (TRUE %in% is.na(overlap)) {
    overlap <- overlap[-is.na(overlap)]
  }
  if (length(overlap)>0) {
    referenceGenome <- referenceGenome[-overlap,]
    assign("referenceGenome", referenceGenome, envir=get(getEnvironment()))
  }
  invisible()
}


#' accessory method for mapping named chromosomes to their pointers in the reference fasta
#'
#' This method will return the numerical index for the named chromosome within the referenceGenomeSequence
#'
#' @param chrId is the identifier for the chromosomes to be dropped
#' @return None
#'
#' @examples
#' \dontrun{
#' getStringSetId("1")
#' }
#'
#' @export
getStringSetId <- function(chrId) {
  if (!(exists("referenceGenome", envir=get(getEnvironment())) & exists("referenceGenomeSequence", envir=get(getEnvironment())))) {
    loadReferenceGenome()
  }
  referenceGenome <- get("referenceGenome", envir=get(getEnvironment()))
  return(referenceGenome[match(as.character(chrId), as.character(referenceGenome[,1])),"sid"])
}




#' returns a DNAStringSet object corresponding to specified chromosome from reference genome
#'
#' This method will return a DNAStringSet object from the reference genome; requires a numeric pointer to the object
#'
#' @param dnaStringSetId is the pointer to use
#' @return DNAStringSet corresponding to pointer provided
#'
#' @examples
#' \dontrun{
#' getChromosomeSequence(1)
#' }
#'
#' @seealso [getStringSetId()] for method to prepare numeric pointer
#' @export
getChromosomeSequence <- function(dnaStringSetId) {
  if (!(exists("referenceGenome", envir=get(getEnvironment())) & exists("referenceGenomeSequence", envir=get(getEnvironment())))) {
    loadReferenceGenome()
  }
  referenceGenomeSequence <- get("referenceGenomeSequence", envir=get(getEnvironment()))
  return(referenceGenomeSequence[[dnaStringSetId]])
}




#' get chromosome identifiers from reference genome
#'
#' This method will return a vector of chromsome identifiers
#'
#' @return vector of names
#'
#' @examples
#' \dontrun{
#' getChromosomeIds("1")
#' }
#'
#' @export
getChromosomeIds <- function() {
  if (!(exists("referenceGenome", envir=get(getEnvironment())) & exists("referenceGenomeSequence", envir=get(getEnvironment())))) {
    loadReferenceGenome()
  }
  return(get("referenceGenome", envir=get(getEnvironment()))[,1])
}



