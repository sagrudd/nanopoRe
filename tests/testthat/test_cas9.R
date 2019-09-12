context("Cas9 and enrichment functionality")

init()
yamlFile <- system.file("extdata", "cas9_demo.yaml", package = "nanopoRe")
referenceGenome <- system.file("extdata", "cas9_demo_ref.fasta", package = "nanopoRe")
bamFile <- system.file("extdata", "cas9_FAK76554.bam", package = "nanopoRe")
unmappedQ <- system.file("extdata", "cas9_FAK76554.unmapped.quals", package = "nanopoRe")
bedTargets <- system.file("extdata", "cas9_demo_target.bed", package = "nanopoRe")

test_that("cas9 - yaml parsing", {
    expect_silent(sourceCas9Parameters(yamlFile))
    expect_silent(yaml <- cas9ParametersToYAML())
    expect_equal(class(yaml), "character")
    expect_silent(html <- cas9ParametersToYAML("Kable"))
    expect_true("kableExtra" %in% class(html))
    expect_null(cas9ParametersToYAML("Excel"))
    expect_vector(getCas9ParameterFields())
    expect_equal(getCachedYAMLValue(field="fastq"), "RawData/FAK76554.fastq.gz")
    expect_null(getCachedYAMLValue(field="sheep")) # nonsense check
})

test_that("cas9 - update config parameters", {
    expect_equal(getCachedYAMLValue(field="reference_genome"), "http://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/cas9_demo_ref.fasta.gz")
    expect_true(setCas9ParameterValue("reference_genome", referenceGenome))
    expect_false(setCas9ParameterValue("nonsensefield", referenceGenome))
    expect_equal(getCachedYAMLValue(field="reference_genome"), referenceGenome)
    expect_false(hasCas9ParameterField("bam_file"))
    expect_true(addCas9ParameterValue("bam_file", bamFile))
    expect_true(hasCas9ParameterField("bam_file"))
    expect_false(addCas9ParameterValue("bam_file", bamFile))
})


test_that("cas9 - enrichment analysis workflow", {
    expect_silent(sourceCas9Parameters(yamlFile))
    expect_false(is_casDataRun())

    expect_true(is_casDataRun())
})


test_that("cas9 - unmapped read q investigation", {
    expect_silent(sourceCas9Parameters(yamlFile))
    addCas9ParameterValue("unmapped_quals", unmappedQ)

    unmapped <- harvestUnmappedQuals(unmappedQ, force=TRUE)
    expect_equal(dim(unmapped), c(99, 2))
})


test_that("parseGenomeCoordinates - test core functionality", {
    expect_silent(sourceCas9Parameters(yamlFile))
    expect_true(setCas9ParameterValue("reference_genome", referenceGenome))
    expect_true(setCas9ParameterValue("target_regions", bedTargets))
    setReferenceGenome(getCachedYAMLValue(field="reference_genome"))
    parseGenomeCoordinates(getCachedYAMLValue(field="target_regions"),
                           getCachedYAMLValue(field="target_proximity"))
    message("Hello World")
})


