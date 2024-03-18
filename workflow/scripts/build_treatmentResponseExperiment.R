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

load("rdata_files/build_treatmentResponseExperiment.RData")

# Should load a variable called doseResponseCurves_TFRI_TNBC_UHN
load(INPUT$treatmentResponse)

treatmentResponse <- doseResponseCurves_TFRI_TNBC_UHN

treatmentNames <- names(treatmentResponse)

raw <- lapply(treatmentNames, function(treatment){
    # message(paste("Processing: ", treatment, sep = "\n\t"))
    samples <- names(treatmentResponse[[treatment]])

    dt <- lapply(samples, function(sample){
        # message(paste("Processing: ", sample, sep = "\n\t"))
        Doses <- unlist(treatmentResponse[[treatment]][[sample]][[1]])
        Viabilities <- unlist(treatmentResponse[[treatment]][[sample]][[2]])

        dt <- data.table::data.table(
            treatmentid = rep(treatment, length(Doses)),
            sampleid = rep(sample, length(Doses)),
            dose = Doses,
            viability = Viabilities
        )
        
    }) |> data.table::rbindlist()
    return(dt)
}) |> data.table::rbindlist()

extrapolated <- lapply(treatmentNames, function(treatment){
    # message(paste("Processing: ", treatment, sep = "\n\t"))
    samples <- names(treatmentResponse[[treatment]])

    dt <- lapply(samples, function(sample){
        # message(paste("Processing: ", sample, sep = "\n\t"))
        Doses <- unlist(treatmentResponse[[treatment]][[sample]][[1]])
        Viabilities <- unlist(treatmentResponse[[treatment]][[sample]][[2]])

        Extrapolated_Doses <- unlist(treatmentResponse[[treatment]][[sample]][[3]][[1]][[1]])
        Extrapolated_Viabilities <- unlist(treatmentResponse[[treatment]][[sample]][[3]][[1]][[2]])

        dt <- data.table::data.table(
            treatmentid = rep(treatment, length(Extrapolated_Doses)),
            sampleid = rep(sample, length(Extrapolated_Doses)),
            Dose = Extrapolated_Doses,
            Viability = Extrapolated_Viabilities
        )
        
    }) |> data.table::rbindlist()
    return(dt)
}) |> data.table::rbindlist()



tdm <- CoreGx::TREDataMapper(rawdata=raw)


CoreGx::rowDataMap(tdm) <- list(
    id_columns = c("treatmentid", "dose"),
    mapped_columns = c()
)

CoreGx::colDataMap(tdm) <- list(
    id_columns = c("sampleid"),
    mapped_columns = c()
)

CoreGx::assayMap(tdm) <- list(
    sensitivity = list(
        id_columns = c("treatmentid", "sampleid", "dose"),
        mapped_columns = c("viability")
    )
)
(uhn_tre <- CoreGx::metaConstruct(tdm))
uhn_tre$sensitivity


uhn_tre_fit <- uhn_tre |> CoreGx::endoaggregate(
    {  # the entire code block is evaluated for each group in our group by
        # 1. fit a log logistic curve over the dose range
        fit <- PharmacoGx::logLogisticRegression(dose, viability,
            viability_as_pct=TRUE)
        # 2. compute curve summary metrics
        ic50 <- PharmacoGx::computeIC50(dose, Hill_fit=fit)
        aac <- PharmacoGx::computeAUC(dose, Hill_fit=fit)
        # 3. assemble the results into a list, each item will become a
        #   column in the target assay.
        list(
            HS=fit[["HS"]],
            E_inf = fit[["E_inf"]],
            EC50 = fit[["EC50"]],
            Rsq=as.numeric(unlist(attributes(fit))),
            aac_recomputed=aac,
            ic50_recomputed=ic50
        )
    },
    assay="sensitivity",
    target="profiles",
    enlist=FALSE,  # this option enables the use of a code block for aggregation
    by=c("treatmentid", "sampleid"),
    nthread=THREADS  # parallelize over multiple cores to speed up the computation
)


uhn_tre_fit$dose_response_curves <- extrapolated



message("Saving the treatment response experiment object to ", OUTPUT$tre)
saveRDS(uhn_tre_fit, file = OUTPUT$tre)
