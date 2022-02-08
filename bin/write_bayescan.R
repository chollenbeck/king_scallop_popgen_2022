write_bscan <- function(genind, file) {
  
  # Make a locus and allele index
  loc_index <- data.frame(locus = adegenet::locNames(genind), loc_code = 1:adegenet::nLoc(genind))
  
  allele_index <- genind@tab %>%
    colnames() %>%
    data.frame(locus.allele = .) %>%
    tidyr::separate(locus.allele, into = c("locus", "allele"), sep = "\\.") %>%
    dplyr::mutate(allele_code = as.integer(allele))
  
  pop_index <- adegenet::popNames(genind) %>%
    data.frame(pop = .) %>%
    dplyr::mutate(pop = gsub(x = pop, pattern = "_\\d+", replacement = "")) %>%
    dplyr::mutate(pop_code = row_number())
  
  # Internal function to get allele counts for each population
  get_bscan_counts <- function(gen, loc_index, allele_index) {
    
    count_tbl <- gen@tab %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "ind") %>%
      tidyr::gather(locus.allele, has.allele, -ind) %>%
      tidyr::separate(locus.allele, into = c("locus", "allele"), sep = "\\.") %>%
      dplyr::group_by(locus, allele) %>%
      dplyr::summarise(count = sum(has.allele, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(loc_index, by = "locus") %>%
      dplyr::left_join(allele_index, by = c("locus", "allele")) %>%
      dplyr::arrange(loc_code, allele_code) %>%
      dplyr::group_by(loc_code) %>%
      dplyr::summarise(n = sum(count), n_alleles = length(count), counts = paste(count, collapse = " "))
    
    return(count_tbl)
    
  }
  
  bscan_counts <- genind %>%
    adegenet::seppop() %>%
    purrr::map(~get_bscan_counts(., loc_index, allele_index))
  
  loc_str <- glue::glue("[loci]=", adegenet::nLoc(genind), "\n\n")
  readr::write_file(x = loc_str, path = file, append = FALSE)
  pop_str <- glue::glue("[populations]=", adegenet::nPop(genind), "\n\n")
  readr::write_file(x = pop_str, path = file, append = TRUE)
  for (i in 1:adegenet::nPop(genind)) {
  
    readr::write_file(x = glue::glue("[pop]={i}\n\n"), path = file, append = TRUE)
    write.table(bscan_counts[[i]],
                file = file,
                quote = FALSE,
                sep = " ",
                col.names = FALSE,
                row.names = FALSE,
                append = TRUE)
    
  }
  
}

