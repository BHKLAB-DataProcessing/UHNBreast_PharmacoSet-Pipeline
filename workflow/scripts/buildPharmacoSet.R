## ------------------- Parse Snakemake Object ------------------- ##
# Check if the "snakemake" object exists
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    
    # setup logger if log file is provided
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)

    save.image(
        file.path("rdata_files/", paste0(snakemake@rule, ".RData"))
    )
}
load("rdata_files/buildPharmacoSet.RData")


tre <- readRDS(INPUT$tre)
mae <- readRDS(INPUT$mae)


# 1.0 Create additional metadata 
# ------------------------------
mae_sampleNames <- lapply(MultiAssayExperiment::colnames(mae), unique) |> 
    unlist() |> 
    unique()

tre_sampleNames <- gsub("(-| |PT)", "", toupper(colnames(tre))) |> 
    unique()

sampleNames <- c(
    lapply(MultiAssayExperiment::colnames(mae), unique) |> 
        unlist() |> 
        unique(),
    colnames(tre) |> 
        unique()
) |> unique()

sampleMetadata <- data.table::fread(INPUT$sampleMetadata)

sampleMetadata$sampleid <- sampleMetadata$cell_line
data.table::setkeyv(sampleMetadata, "sampleid")
sampleMetadata <- sampleMetadata[sampleid %in% sampleNames,] |> unique()


# subset the tre for the ones we have sample info for
tre <- tre[, sampleMetadata$sampleid]

sample <- as.data.frame(
    sampleMetadata, 
    row.names = sampleMetadata$sampleid
)


treatmentNames <- rowData(tre)[, unique(treatmentid)]

treatment <- data.frame(
    treatmentid = treatmentNames,
    row.names = treatmentNames
)

name <- "UHNBreast"


pset <- PharmacoGx::PharmacoSet2(
    name = name,
    treatment = treatment,
    sample = sample,
    molecularProfiles = mae,
    treatmentResponse = tre,
    perturbation = list(),
    curation = list(
        sample = sample, 
        treatment = treatment, 
        tissue = data.frame()),
    datasetType = "sensitivity"
)

message(paste(capture.output(show(pset)), collapse = "\n\t"))

message("Object Size (pset):")
object.size(pset) |> print(unit = "auto")
print(paste("Saving PharmacoSet object to", OUTPUT[[1]]))
dir.create(dirname(OUTPUT[[1]]), recursive = TRUE, showWarnings = FALSE)
saveRDS(pset, file = OUTPUT[[1]])
