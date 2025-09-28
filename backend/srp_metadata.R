library(rentrez)
library(xml2)
library(jsonlite)

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
        
        study_info <- list(
          accession = accession,
          study_accession = xml_attr(xml_find_first(doc, "//STUDY"), "accession"),
          title = safe_extract_xml2("//STUDY_TITLE", doc),
          abstract = safe_extract_xml2("//STUDY_ABSTRACT", doc),
          description = safe_extract_xml2("//STUDY_DESCRIPTION", doc),
          study_type = xml_attr(xml_find_first(doc, "//STUDY_TYPE"), "existing_study_type"),
          center = xml_attr(xml_find_first(doc, "//STUDY"), "center_name"),
          submission_date = xml_attr(xml_find_first(doc, "//STUDY"), "created")
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
      study_accession = ifelse(is.null(study$study_accession) || is.na(study$study_accession), 
                              acc, study$study_accession),
      title = ifelse(is.null(study$title), NA, study$title),
      abstract = ifelse(is.null(study$abstract), NA, study$abstract),
      description = ifelse(is.null(study$description), NA, study$description),
      study_type = ifelse(is.null(study$study_type), NA, study$study_type),
      center = ifelse(is.null(study$center), NA, study$center),
      submission_date = ifelse(is.null(study$submission_date), NA, study$submission_date),
      error = ifelse(is.null(study$error), NA, study$error),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, df_list)
}

harvest_bundle_metadata <- function(org) {
  bundle_path <- paste("/mnt/md0/dee2/sradb/big_proj/",org,sep="")
  zips <- list.files(bundle_path,pattern="zip$")
  srps <- sapply(strsplit(zips,"_"),"[[",1)
  zips <- zips[which(! duplicated(srps))] # exclude duplicate bundles without GSE info
  srps <- srps[which(! duplicated(srps))]
  gse <- gsub(".zip","",sapply(strsplit(zips,"_"),"[[",2))
  myfilename <- paste("../sradb/",org,"_srp.tsv",sep="")
  file.copy(myfilename,"tmp0",overwrite=TRUE)
  system("cut -f1 tmp0  > tmp1")
  srp_existing <- gsub('"',"",readLines("tmp1"))
  srps_new <- setdiff(srps,srp_existing)
  res2 <- get_sra_studies_metadata_xml2(srps_new)
  res2df <- metadata_to_dataframe(res2)
  res2df$GSE <- gse[which(srps %in% res2df$query_accession)]
  zips2 <- zips[which(srps %in% res2df$query_accession)] # only include bundles with metadata
  res2df$URL <- paste("https://dee2.io/huge/",org,"/",zips2,sep="")
  write.table(x=res2df,file=myfilename,sep="\t",row.names=FALSE, append=TRUE)
  SYSCOMMAND <- gsub("myorg",org,"scp -i ~/.ssh/dee2_2025 ../sradb/myorg_srp.tsv ubuntu@dee2.io:/dee2_data/metadata/")
  system(SYSCOMMAND)
  unlink("tmp1")
  unlink("tmp0")
}

#orgs <- c("athaliana", "bdistachyon", "celegans", "dmelanogaster", "drerio",
#  "ecoli", "gmax", "hsapiens", "hvulgare", "mmusculus", "osativa", "ptrichocarpa",
#  "rnorvegicus", "sbicolor", "scerevisiae", "slycopersicum", "stuberosum",
#  "taestivum", "vvinifera", "zmays")
#lapply(orgs,harvest_bundle_metadata)
harvest_bundle_metadata(org)
