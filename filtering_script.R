library(tidyverse)
library(DESeq2)
library(ggplot2)
print("Filtering script loaded")
# Base path for saving results
base_output_path <- "C:/Users/Aditi/Documents/Desktop/R/DGE Analysis Main/results"
# Create output directory if it doesn't exist
if (!dir.exists(base_output_path)) {
  dir.create(base_output_path, recursive = TRUE)
}

###FUNCTION 1
create_deseq_formula <- function(meta) {
  # Convert all columns to factors
  for (col in names(meta)) {
    meta[[col]] <- factor(meta[[col]])
  }
  
  priority_columns <- c(
    "Timepoint..after.stressor.", 
    "Replicate" 
  )
  
  # Conditionally add Cell.Line if it has more than one unique value
  if (length(unique(meta$Cell.Line)) > 1) {
    priority_columns <- c(priority_columns, "Cell.Line")
  }
  
  # Select only the priority columns that exist in the metadata
  selected_columns <- intersect(priority_columns, names(meta))
  
  # Fallback to first column if no priority columns found
  if (length(selected_columns) == 0) {
    selected_columns <- names(meta)[1]
  }
  
  # Create formula dynamically
  formula <- as.formula(paste("~", paste(selected_columns, collapse = " + ")))
  
  return(formula)
}

process_gse <- function(gse_id, dataset_meta) {
  print(paste("Processing GEO ID:", gse_id))
  
  # Create a folder for this specific GSE ID if it doesn't exist
  gse_output_path <- file.path(base_output_path, gse_id)
  if (!dir.exists(gse_output_path)) {
    dir.create(gse_output_path, recursive = TRUE)
  }
  
  # Step 1: Load Metadata
  meta <- dataset_meta %>% 
    filter(GSE.ID == gse_id) %>% 
    tibble::column_to_rownames(var = "Sample.ID")
  
  # Step 2: Load Raw Counts
  raw_counts_path <- file.path("raw_counts", paste0(gse_id, "_raw_counts_GRCh38.p13_NCBI.tsv.tsv"))
  if (!file.exists(raw_counts_path)) stop(paste("Raw counts file not found:", raw_counts_path))
  raw_counts <- read.table(raw_counts_path, header = TRUE, row.names = 1)
  
  # Step 3: Filter Raw Counts to Match Metadata
  matching_samples <- intersect(colnames(raw_counts), rownames(meta))
  counts_filtered <- raw_counts[, matching_samples]
  meta_filtered <- meta[matching_samples, , drop = FALSE]
  
  # Step 4: Perform CPM Filtering
  total_counts <- colSums(counts_filtered)
  cpm_data <- sweep(counts_filtered, 2, total_counts, "/") * 1e6
  keep <- rowSums(cpm_data > 1) >= (ncol(counts_filtered) * 0.5)
  counts_filtered <- counts_filtered[keep, ]
  
  # Step 5: Dynamic DESeq2 formula creation
  deseq_formula <- create_deseq_formula(meta_filtered)
  
  if (all(dataset_meta$single.timecourse[dataset_meta$GSE.ID == gse_id] == "single")) {
    # Design for single condition
    dds <- DESeqDataSetFromMatrix(
      countData = counts_filtered, 
      colData = meta_filtered, 
      design = ~ Replicate + Stressor..with.concentration.
    )
    
    # Find unique Control label based on Control.Treatment
    control_label <- unique(meta_filtered$Stressor..with.concentration.[meta_filtered$Control.Treatment == "Control"])
    
    if (length(control_label) != 1) stop("Ambiguous 'Control' label found. Ensure one unique control label exists.")
    
    # Relevel to make Control the reference
    dds$Stressor..with.concentration. <- relevel(factor(dds$Stressor..with.concentration.), ref = control_label)
    
    # Get unique treatments excluding the control
    unique_treatments <- unique(meta_filtered$Stressor..with.concentration[meta_filtered$Stressor..with.concentration != control_label])
    
    # Run DESeq2 contrasts for each treatment
    contrast_results <- lapply(unique_treatments, function(treatment) {
      results(dds, contrast = c("Stressor..with.concentration.", treatment, control_label))
    })
    
    # Assign names to the contrast results
    names(contrast_results) <- unique_treatments
    
    # Save each contrast result individually in the GSE-specific folder
    lapply(names(contrast_results), function(treatment) {
      # Create a file path for each contrast result
      file_path <- file.path(gse_output_path, paste0(gse_id, "_", treatment, "_contrast.RData"))
      
      # Save each contrast result as a separate RData file
      save(contrast_results[[treatment]], file = file_path)
    })
    
  } 
  
  if (all(dataset_meta$single.timecourse[dataset_meta$GSE.ID == gse_id] == "timecourse")) {
    # Convert all columns to factors
    for (col in names(meta_filtered)){
      meta_filtered[[col]] <- factor(meta_filtered[[col]])
    }
    
    # Define a function to dynamically reduce the design formula
    reduce_design_formula <- function(deseq_formula, base_variable = "Timepoint..after.stressor.") {
      # Extract the terms in the design formula (remove intercept ~1)
      terms <- all.vars(deseq_formula)[-1]  # Removing the intercept term
      
      # If only Timepoint is present, return formula ~ 1 (intercept-only)
      if (length(terms) == 1 && terms == base_variable) {
        reduced_formula <- as.formula("~ 1")
      } else {
        # Remove the base variable (Timepoint..after.stressor.)
        reduced_formula <- update(deseq_formula, paste("~ . -", base_variable))
      }
      
      return(reduced_formula)
    }
    
    # Call the function to reduce the formula dynamically
    reduced_formula <- reduce_design_formula(deseq_formula)
    
    # Apply the reduced formula in the DESeq2 pipeline
    dds <- DESeqDataSetFromMatrix(
      countData = counts_filtered, 
      colData = meta_filtered, 
      design = deseq_formula
    )
    
    # Perform Likelihood Ratio Test
    dds_lrt <- DESeq(dds, test = "LRT", reduced = reduced_formula)
    res_LRT <- results(dds_lrt)
    
    # Filter significant genes
    res_LRT_tb <- res_LRT %>%
      data.frame() %>%
      rownames_to_column(var = "gene") %>%
      as_tibble()
    
    sigLRT_genes <- res_LRT_tb %>%
      filter(padj < 0.01, baseMean > 10, abs(log2FoldChange) > 1)
    
    # Save significant LRT genes in the GSE-specific folder
    save(sigLRT_genes, file = file.path(gse_output_path, paste0(gse_id, "_LRT_genes.RData")))
  }
  # Step 7: Normalization
  dds <- estimateSizeFactors(dds)
  normalized_counts <- counts(dds, normalized = TRUE)
  vsd <- vst(dds, blind = TRUE)
  
  # Step 8: PCA Plot
  # Determine intgroup based on specific columns
  # Check for Timepoint column (with the specific name you mentioned)
  intgroup <- c()
  
  if ("Timepoint..after.stressor." %in% names(meta_filtered)) {
    intgroup <- c(intgroup, "Timepoint..after.stressor.")
  }
  
  # Check for Cell Line column
  if ("Cell.Line" %in% names(meta_filtered)) {
     intgroup <- c(intgroup, "Cell.Line")
  }
  # If no columns found, use the first categorical column as a fallback
  if (length(intgroup) == 0) {
    intgroup <- names(meta_filtered)[sapply(meta_filtered, is.factor)][1]
  }
  
  # Create PCA plot using the determined intgroup
  tryCatch({
    pca_plot <- plotPCA(vsd, intgroup = intgroup) +
      labs(title = paste("PCA Plot -", gse_id)) +
      theme_minimal()
  }, error = function(e) {
    warning(paste("Could not create PCA plot for", gse_id, ":", e$message))
    pca_plot <- NULL
  })
  
  # Create output directory
  gse_output_dir <- file.path("results", gse_id)
  dir.create(gse_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save results
  write.table(
    counts_filtered, 
    file.path(gse_output_dir, paste0(gse_id, "_filtered_counts.tsv")), 
    sep = "\t", 
    quote = FALSE
  )
  
  write.table(
    normalized_counts, 
    file.path(gse_output_dir, paste0(gse_id, "_normalized_counts.tsv")), 
    sep = "\t", 
    quote = FALSE
  )
  
  # Save the DESeq2 object
  save(dds, file = file.path(gse_output_dir, paste0(gse_id, "_dds_object.RData")))
  
  # Save PCA plot if created
  if (!is.null(pca_plot)) {
    ggsave(
      filename = file.path(gse_output_dir, paste0(gse_id, "_pca_plot.png")), 
      plot = pca_plot, 
      width = 10, 
      height = 8
    )
  # Save DESeq2 object in the GSE-specific folder
  save(dds, file = file.path(gse_output_path, paste0(gse_id, "_dds_object.RData")))
  
  #Save filtered counts
  save(counts_filtered, file = file.path(gse_output_path, paste0(gse_id, "_filtered_counts.RData")))
  
  return(list(
    raw_counts = counts_filtered,
    metadata = meta_filtered,
    deseq_object = dds,
    normalized_counts = normalized_counts,
    pca_plot = pca_plot
  ))
}
  }
