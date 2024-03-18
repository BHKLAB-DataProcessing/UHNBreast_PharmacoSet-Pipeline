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
load("rdata_files/build_MultiAssayExperiment.RData")

sampleMetadata <- data.table::fread(INPUT$sampleMetadata)

# this should load a variable called "rse_list"
load(INPUT$rse_list)


colData <- as.data.frame(sampleMetadata, row.names = sampleMetadata$Run)

colData$sampleid <- colData$Run
colData$batchid <- 1



sampleNames <- lapply(
    rse_list, SummarizedExperiment::colnames
    ) |> 
    unlist() |> 
    unique()




if(!all(sampleNames %in% sampleMetadata$Run)){
    print("Not all sample names are in the sample metadata")

    message("Removing the following sample names from the summarized experiments:\n\t")
    missing <- setdiff(se_sampleNames, sampleMetadata$Run)
    message("\t",paste(missing, collapse = "\n\t"))

    se_list <- lapply(se_list, function(x){
        x[,!colnames(x) %in% missing]
    })
}

summarizedExperimentLists <- sapply(rse_list, function(x){
    x@colData <- MultiAssayExperiment::DataFrame(
        sampleid = colnames(x),
        batchid = rep(NA, ncol(x)),
        row.names = colnames(x)
    )
    x
})
ExpList <- MultiAssayExperiment::ExperimentList(summarizedExperimentLists)
message(paste("ExperimentList:", capture.output(show(ExpList)), sep = "\n\t"))


# Create a sample map for each experiment in the ExperimentList
sampleMapList <- lapply(summarizedExperimentLists, function(se){
    data.frame(
        primary = colnames(se),
        colname = colnames(se),
        stringsAsFactors = FALSE
    )
})
names(sampleMapList) <- names(ExpList)
message(paste("Sample map list:", capture.output(str(sampleMapList)), sep = "\n\t"))


mae <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = ExpList,
    colData = colData,
    sampleMap = MultiAssayExperiment::listToMap(sampleMapList)
)

# Write Output
# ------------
message("Saving MultiAssayExperiment to: ", OUTPUT$mae)
saveRDS(mae, file = OUTPUT$mae)
