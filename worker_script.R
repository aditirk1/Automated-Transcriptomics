library(httr)
library(R.utils)
# Function to download and extract raw count files
download_and_extract_counts <- function(geo_ids, output_dir = "raw_counts") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (geo_id in geo_ids) {
    # Construct the full download URL for raw counts
    full_url <- paste0(
      "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=", 
      geo_id, 
      "&format=file&file=", 
      geo_id, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    )
    
    # Destination file path
    gz_file <- file.path(output_dir, paste0(geo_id, "_raw_counts_GRCh38.p13_NCBI.tsv.gz"))
    
    # Attempt download
    tryCatch({
      # Download the file
      GET(full_url, write_disk(gz_file, overwrite = TRUE))
      cat("Successfully downloaded raw counts for", geo_id, "\n")
      
      # Extract the .gz file to .tsv format
      tsv_file <- sub(".gz$", ".tsv", gz_file)
      gunzip(gz_file, destname = tsv_file, remove = TRUE)
      cat("Successfully extracted the file:", tsv_file, "\n")
      
    }, error = function(e) {
      warning("Error downloading or extracting file for ", geo_id, ": ", e$message)
    })
  }
}