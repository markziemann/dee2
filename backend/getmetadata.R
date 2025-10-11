library(rentrez)
library(xml2)
library(jsonlite)
library(kableExtra)

#Functions
get_sra_studies_metadata_xml2 <- function(study_accessions) {
  results <- list()
  
  for (accession in study_accessions) {
    cat("Fetching:", accession, "\n")
    
    tryCatch({
      # Search for the study
      search_result <- entrez_search(
        db = "sra",
        term = paste0(accession, "[Accession]"),
        retmax = 100
      )
      
      if (length(search_result$ids) > 0) {
        # Fetch XML data
        xml_data <- entrez_fetch(
          db = "sra",
          id = search_result$ids[1],
          rettype = "xml"
        )
        
        # Parse XML with xml2
        doc <- read_xml(xml_data)
        
        # Helper function to safely extract values
        safe_extract_xml2 <- function(xpath_expr, doc) {
          result <- xml_find_all(doc, xpath_expr)
          if (length(result) == 0) {
            return(NA)
          } else if (grepl("@", xpath_expr)) {
            # It's an attribute
            return(xml_text(result[1]))
          } else {
            # It's an element
            return(xml_text(result[1]))
          }
        }

        get_date <- function(doc) {
          # Try various date fields in order of preference
          dates <- list(
            study_published = xml_attr(xml_find_first(doc, "//STUDY"), "published"),
            study_created = xml_attr(xml_find_first(doc, "//STUDY"), "created"),
            study_updated = xml_attr(xml_find_first(doc, "//STUDY"), "updated"),
            submission_published = xml_attr(xml_find_first(doc, "//SUBMISSION"), "published"),
            submission_created = xml_attr(xml_find_first(doc, "//SUBMISSION"), "created"),
            package_published = xml_attr(xml_find_first(doc, "//EXPERIMENT_PACKAGE"), "published"),
            run_published = xml_attr(xml_find_first(doc, "//RUN"), "published")
          )
          
          # Return the first non-NA date
          for (date_name in names(dates)) {
            if (!is.na(dates[[date_name]])) {
              return(list(date = dates[[date_name]], source = date_name))
            }
          }
          return(list(date = NA, source = NA))
        }

        date_info <- get_date(doc)

        study_info <- list(
          accession = accession,
          study_accession = xml_attr(xml_find_first(doc, "//STUDY"), "accession"),          
          species = safe_extract_xml2("//SCIENTIFIC_NAME", doc),
          title = safe_extract_xml2("//STUDY_TITLE", doc),
          study_type = xml_attr(xml_find_first(doc, "//STUDY_TYPE"), "existing_study_type"),
          date_published = sapply(strsplit(date_info[[1]]," "),"[[",1)
        )
     
        results[[accession]] <- study_info
        
      } else {
        cat("No results found for", accession, "\n")
        results[[accession]] <- list(accession = accession, error = "No results found")
      }
      
    }, error = function(e) {
      cat("Error processing", accession, ":", e$message, "\n")
      results[[accession]] <- list(accession = accession, error = e$message)
    })
    
    Sys.sleep(3)  # Rate limiting
  }
  
  return(results)
}

metadata_to_dataframe <- function(metadata_list) {
  df_list <- lapply(names(metadata_list), function(acc) {
    study <- metadata_list[[acc]]
    data.frame(
      query_accession = acc,
      species = ifelse(is.null(study$species), NA, study$species),
      title = ifelse(is.null(study$title), NA, study$title),
      study_type = ifelse(is.null(study$study_type), NA, study$study_type),
      submission_date = ifelse(is.null(study$date_published), NA, study$date_published),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, df_list)
}


# work
myfiles <- list.files(".",pattern="confirmed")
srps <- sapply(strsplit(myfiles,"\\."),"[[",1)
res2 <- get_sra_studies_metadata_xml2(srps)
res2df <- metadata_to_dataframe(res2)
colnames(res2df) <- c("SRA Project Accession","Species","Project Title","Study Type","Date Published")
res2df$`Date Processed` <- unlist(lapply(myfiles,function(x) {
  ctime <- file.info(x)$ctime
  sapply(strsplit(as.character(ctime)," "),"[[",1) }
) )
res2df <- res2df[order(res2df$`Date Processed`,decreasing=TRUE),]
saveRDS(res2df,"metadata.rds")
writeLines(kbl(res2df,row.names=FALSE),con="tbl.html")
writeLines(kbl(head(res2df,5),row.names=FALSE),con="top.html")
system("scp -i ~/.ssh/dee2_2025 *html ubuntu@dee2.io:~/dee2/frontend/html/")

