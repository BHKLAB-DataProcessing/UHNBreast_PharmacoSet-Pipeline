#####################################################################
####################### Load the environment ########################
######################################################################
# This checks if the snakemake object exists. 
# an alternative if running as a script
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    if(length(snakemake@log)>0) 
        sink(snakemake@log[[1]], FALSE, c("output", "message"), TRUE)
}else{
    INPUT <- list(
        annotation_file = "references/GRCh38_v45/annotation.gtf",
        quant_files = c(
            "procdata/rnaseq/kallisto_v0.46.1_GRCh38.45/SRR2532365/abundance.h5",
            "procdata/rnaseq/kallisto_v0.46.1_GRCh38.45/SRR2532333/abundance.h5"
        )
    )
    OUTPUT <- list(
        rse_list = "procdata/rnaseq/kallisto_v0.46.1_GRCh38.45/rse_list.RData"
    )
    WILDCARDS <- list(
        tool = "kallisto",
        tool_version = "0.46.1",
        ref_build = "GRCh38",
        ref_version = "45"
    )
    THREADS <- 1
}
#####################################################################
####################### Generalized functions #######################
#####################################################################

#' getTranscripts Function
#'
#' This function imports transcript-level quantification files and performs necessary preprocessing steps.
#' Written to handle the transcription annotation files from GENCODE and clean the rownames. 
#'
#' @param quant_files A character vector specifying the paths to the quantification files.
#' @param tx2gene A TxDb object or a character vector specifying the path to the transcript-to-gene mapping file.
#' @param tool A character string specifying the quantification tool used (e.g., "salmon", "kallisto").
#' @param countsFromAbundance A character string specifying whether to derive counts from abundance values ("yes" or "no").
#'
getTranscripts <- function(quant_files, tx2gene, tool, countsFromAbundance = "no"){
    transcripts <- tximport::tximport(
        files = quant_files, 
        type = tool, 
        txOut = TRUE, 
        ignoreTxVersion = FALSE,
        countsFromAbundance=countsFromAbundance
    )
    rownames(transcripts$counts) <- sub("\\|.*", "", rownames(transcripts$counts))
    rownames(transcripts$abundance) <- sub("\\|.*", "", rownames(transcripts$abundance))
    rownames(transcripts$length) <- sub("\\|.*", "", rownames(transcripts$length))

    return(transcripts)
}

#' getGenes Function
#'
#' This function imports gene expression quantification files and performs gene-level summarization using the tximport package.
#'
#' @param quant_files A character vector specifying the paths to the quantification files.
#' @param tx2gene A TxDb object or a named character vector mapping transcript IDs to gene IDs.
#' @param tool A character string specifying the quantification tool used (e.g., "salmon", "kallisto").
#'
getGenes <- function(quant_files, tx2gene, tool){
    print(sprintf("Loading %s Gene Data for %s number of files", tool, length(quant_files)))
    genes <- tximport::tximport(
        quant_files, 
        type = tool, 
        tx2gene = tx2gene, 
        txIn = ifelse(tool == "rsem", FALSE, TRUE), 
        ignoreAfterBar = TRUE, 
        ignoreTxVersion = FALSE
    )
    return(genes)
}


####################################################################
####################### Load libraries #############################
####################################################################

# set lib_path to "/usr/local/lib/R/site-library"
# TODO:: find a way to circumvent the hardcoding of the lib_path
.libPaths(c("/usr/local/lib/R/site-library", "/cluster/home/t119797uhn/R/x86_64-pc-linux-gnu-library/4.2"))
library(tximport)
library(GenomicFeatures)
library(rtracklayer)
library(AnnotationDbi)


#####################################################################
####################### Load Data and Process #######################
#####################################################################
# 1. Load the genomic annotation

print("Loading Genomic annotation Data")
txdb <- GenomicFeatures::makeTxDbFromGFF(file = INPUT$annotation_file)
keys <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, keys, "GENEID", "TXNAME")
GRanges <- rtracklayer::import(INPUT$annotation_file)


# 2. Load the quantification files
# Tool should be one of "salmon", "kallisto" or "rsem"
tool <- WILDCARDS$tool
stopifnot(tool %in% c("salmon", "kallisto", "rsem"))

quant_files <- INPUT$quant_files
names(quant_files) <- basename(dirname(quant_files))

if(tool %in% c("salmon", "kallisto")){
    rnaseq.genes<- getGenes(quant_files, tx2gene, tool)
    rnaseq.transcripts <- getTranscripts(quant_files, tx2gene, tool, countsFromAbundance = "no")
}
### Save the data
assays <- list(
    genes = rnaseq.genes$abundance,
    genes_counts = rnaseq.genes$counts,
    genes_length = rnaseq.genes$length,
    transcripts = rnaseq.transcripts$abundance,
    transcripts_counts = rnaseq.transcripts$counts,
    transcripts_length = rnaseq.transcripts$length
)
str(assays)

metadata <- WILDCARDS[which(names(WILDCARDS) != "")]

print(names(assays))
rse_list <- lapply(names(assays), function(x){
    assay <- assays[[x]]

    if(grepl("genes", x)){
        granges <-  GRanges[GRanges$gene_id %in% rownames(assay) & GRanges$type == "gene",]
        assay <- assay[rownames(assay) %in% granges$gene_id,]

        granges <- granges[match(rownames(assay), granges$gene_id),]
    } else {
        granges <-  GRanges[GRanges$transcript_id %in% rownames(assay) & GRanges$type == "transcript",]
        assay <- assay[rownames(assay) %in% granges$transcript_id,]

        granges <- granges[match(rownames(assay), granges$transcript_id),]
    }

    # Align the rownames of the assay to the rownames of the granges
    

    SummarizedExperiment::SummarizedExperiment(
        assays = list(
            expr = assay
        ),
        rowRanges = granges,
        metadata = metadata
    )
})
names(rse_list) <- names(assays)
print(rse_list)


print(paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " Saving RData to ", OUTPUT$rse_list))
object.size(rse_list)  |> print()

dir.create(dirname(OUTPUT$rse_list), recursive = TRUE, showWarnings = FALSE)
save(rse_list, file = OUTPUT$rse_list, compress = "bzip2", compression_level = 9)
