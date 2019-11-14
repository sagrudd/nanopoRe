
#' Initialise the nanopoRe environment
#'
#' Creates a nanopoRe environment; package specific parameters and values will
#' be stored within this environment; the name of the environment is defined
#' internally
#'
#' @return None
#'
#' @examples
#' init()
#'
#' @export
init <- function() {
    eval(parse(text = paste0(
        nanopoRe.env.name, " <<- new.env(parent=emptyenv())")))
    setRpath(file.path("Analysis", "R"))
}

#' check nanopoRe environment
#'
#' performs a sanity check to ensure that nanopoRe environment is initialised
#'
#' @return a logical defining whether the environment is initialised
isInitialised <- function() {
    eval(parse(text = paste0(
        "return(exists(\"", nanopoRe.env.name, "\", mode=\"environment\"))")))
}
# nanopoRe.env.name <- 'nanopoRe.env' usethis::use_data(
#     nanopoRe.env.name, internal=TRUE)
getEnvironment <- function() {
    if (!isInitialised()) {
        init()
    }
    return(nanopoRe.env.name)
}

getCachedFileObject <- function(objectName, fileName, subenv=NULL) {
    if (hasCachedObject(objectName, subenv)) {
        return(getCachedObject(objectName, subenv))
    } else {
        setCachedObject(objectName, readRDS(file = fileName), subenv)
        return(getCachedObject(objectName, subenv))
    }
}


selectEnvironment <- function(subenv=NULL) {
    envir <- get(getEnvironment())
    if (!is.null(subenv)) {
        envir <- getCachedObject(subenv)
    }
    return(envir)
}


setCachedObject <- function(objectName, data, subenv=NULL) {
    assign(objectName, data, envir=selectEnvironment(subenv))
}

hasCachedObject <- function(objectName, subenv=NULL) {
    if (exists(objectName, envir=selectEnvironment(subenv))) {
        return(TRUE)
    }
    return(FALSE)
}

getCachedObject <- function(objectName, subenv=NULL) {
    return(get(objectName, envir=selectEnvironment(subenv)))

}


#' list objects that are stored within the nanopoRe working space
#'
#' The nanopoRe package makes use of R environments to store data in memory,
#' but within an internally-named environment; this method dumps out (for
#' reasons of QC and transparency) a complete listing of these stored objects
#'
#' @importFrom utils ls.str
#' @param subenv the nested environment to consider
#' @return character vector of objects stored in the nanopoRe environment
#'
#' @examples
#' cachedObjects <- listCachedObjects()
#'
#' @export
listCachedObjects <- function(subenv=NULL) {
    return(ls.str(envir=selectEnvironment(subenv)))
}



#' list the parameter fields stored in a cached YAML config file
#'
#' list the parameter fields stored in a cached YAML config file
#'
#' @param yaml the name of the YAML object
#' @param subenv the sub-environment to process
#' @return vector of field names
#'
#' @examples
#' init()
#' yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
#' importConfigYAML(yamlFile=yamlFile)
#' getCachedYAMLFields()
#'
#' @export
getCachedYAMLFields <- function(yaml="config", subenv=NULL) {
    return(names(getCachedObject(yaml, subenv)))
}


#' check for available named cas9 parameter field
#'
#' check for available named cas9 parameter field
#'
#' @param yaml the name of the YAML object
#' @param field is the name of a parameter field
#' @param subenv the sub-environment to process
#' @return TRUE or FALSE
#'
#' @examples
#' init()
#' yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
#' importConfigYAML(yamlFile=yamlFile)
#' hasCachedYAMLField(field="bam_file")
#' hasCachedYAMLField(field="study_name")
#'
#' @export
hasCachedYAMLField <- function(yaml="config", field, subenv=NULL) {
    params <- getCachedYAMLFields(yaml, subenv)
    return(field %in% params)
}



#' get the stored parameter value from named field
#'
#' get the stored parameter value from named field
#'
#' @param yaml the name of the YAML object
#' @param field the parameter field to lookup
#' @param subenv the sub-environment to process
#' @return value from config
#'
#' @examples
#' init()
#' yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
#' importConfigYAML(yamlFile=yamlFile)
#' getCachedYAMLValue(field="fastq")
#'
#' @export
getCachedYAMLValue <- function(yaml="config", field, subenv=NULL) {
    if (hasCachedYAMLField(yaml, field, subenv)) {
        params <- getCachedObject(yaml, subenv)
        return(params[[field]])
    }
    return(NULL)
}



#' update the stored parameter value from named config field
#'
#' update the stored parameter value from named config field
#'
#' @param yaml the name of the YAML object
#' @param field the parameter field to lookup
#' @param value to value to update
#' @param subenv the sub-environment to process
#' @return TRUE or FALSE if the change can be made
#'
#' @examples
#' init()
#' yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
#' importConfigYAML(yamlFile=yamlFile)
#' referenceGenome <- system.file("extdata", "cas9_demo_ref.fasta",
#'     package = "nanopoRe")
#' setCachedYAMLValue(field="reference_genome", value=referenceGenome)
#' getCachedYAMLValue(field="reference_genome")
#'
#' @export
setCachedYAMLValue <- function(yaml="config", field, value, subenv=NULL) {
    if (hasCachedYAMLField(yaml, field, subenv=subenv)) {
        params <- getCachedObject(yaml, subenv)
        params[[field]] <- value
        setCachedObject(yaml, params, subenv)
        return(TRUE)
    }
    return(FALSE)
}


#' add a stored parameter value for the named config field
#'
#' add a stored parameter value for the named config field
#'
#' @param yaml the name of the yaml container to create
#' @param field the parameter field to lookup
#' @param value to value to update
#' @param subenv the name of subenvironment for holding the yaml
#' @return TRUE or FALSE if the change can be made
#'
#' @examples
#' init()
#' yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
#' bamFile <- system.file("extdata", "cas9_FAK76554.bam", package = "nanopoRe")
#' importConfigYAML(yamlFile=yamlFile)
#' addCachedYAMLValue(field="bam_file", value=bamFile)
#' getCachedYAMLValue(field="bam_file")
#'
#' @export
addCachedYAMLValue <- function(yaml="config", field, value, subenv=NULL) {
    if (!hasCachedYAMLField(yaml, field, subenv)) {
        params <- getCachedObject(yaml, subenv)
        params[[field]] <- value
        setCachedObject(yaml, params, subenv)
        return(TRUE)
    }
    return(FALSE)
}




#' export the stored config parameters parsed from YAML file
#'
#' export the stored config parameters parsed from YAML file
#'
#' @import yaml
#' @param yaml the name of the yaml container to create
#' @param format is the format to return results as (YAML|Kable|list)
#' @param subenv the name of subenvironment for holding the yaml
#' @return YAML object
#'
#' @examples
#' init()
#' yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
#' importConfigYAML(yamlFile=yamlFile)
#' cachedYAMLToYAML(format="markdown")
#'
#' @export
cachedYAMLToYAML <- function(yaml="config", format=NA, subenv=NULL) {
    params <- getCachedObject(yaml, subenv)

    if (is.na(format)) {
        return(as.yaml(params))
    } else {
        table <- kable(t(as.data.frame(params, row.names=NULL)),
            format=format, caption="Configuration parameters",
            booktabs=TRUE, table.envir='table*', linesep="", escape=FALSE) %>%
            kable_styling(c("striped"))
        return(table)
    } else if (format=="list") {
        return(params)
    }
    return(NULL)
}


#' import the cas9 parameters from YAML file
#'
#' import the cas9 parameters from YAML file
#'
#' @import yaml
#' @param yaml the name of the yaml container to create
#' @param yamlFile is a path to a stored YAML config file
#' @param subenv the name of subenvironment for holding the yaml
#' @return None
#'
#' @examples
#' init()
#' yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
#' importConfigYAML(yamlFile=yamlFile)
#'
#' @export
importConfigYAML <- function(yaml="config", yamlFile, subenv=NULL) {
    yamlData <- yaml.load_file(yamlFile)
    if (is.null(subenv)) {
        setCachedObject(yaml, yamlData)
    } else {
        setCachedObject(subenv, new.env())
        setCachedObject(yaml, yamlData, subenv)
    }
}


