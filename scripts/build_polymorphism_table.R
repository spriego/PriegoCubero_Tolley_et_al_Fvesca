#' Build a polymorphism table where each row is a detected polymorphism, and each column ("present", "absent") lists comma-separated accessions with and without the polymorphism, respectively. Also computes the minor allele frequency (MAF).
#'
#' @param polymorphism_file   Path to deletions BED (with columns: chr,start,end,strand,TE_id,accessions_with_polymorphism) or insertions BED (with columns: chr,start,end,chr_te,start_te,end_te,TE_id,present,accessions_with_polymorphism)
#' @param insertions Logical. TRUE for insertions input; FALSE for deletions.
#' @param accessions_all   Character vector of all accessions
#' @param total_accessions Integer, denominator for % missing (default 210)
#' @return A tibble with polymorphism coordinates, list of comma-separated accessions with present or absent polymorphism, and MAF
build_polymorphism_table <- function(
    polymorphism_file,
    insertions = FALSE,
    accessions_all,
    total_accessions = 210
) {
  df <- 
    if (isTRUE(insertions)) {
      read_tsv(
        file = polymorphism_file,
        col_names = c("chr", "start", "end", "chr_te", "start_te", "end_te", "TE_id", "present")
      ) %>%
        select(-chr_te, -start_te, -end_te) %>%
        unite(col = "TE_id", sep = "-", chr:TE_id, remove = FALSE) %>%
        separate(TE_id, sep = ",", into = c("TE_id"), remove = TRUE) %>%
        group_by(TE_id, chr, start, end) %>%
        # Some entries indicate insertions for multiple TEs, separated by commas. 
        # As a result, the refinement step generates two rows 
        # representing insertions for different accessions. Although these insertions 
        # occur at the same locus, the start and end positions differ depending on 
        # whether they refer to the first or the second TE.
        # 
        # To simplify this, I collapse the accessions carrying the polymorphism 
        # and assume that it occurs only in the first TE listed.
        summarise(present = paste(unique(present), collapse = ","), .groups = "drop") %>% 
        ungroup() %>%
        select(-chr, -start, -end)
    }   else {
      read_tsv(
        polymorphism_file,
        col_names = c("chr", "start", "end", "strand", "TE_id", "present")
      )
    }
  
  intermediate_df <- 
    df %>%
    # Add column of accessions not carrying the polymorphism
    mutate(
      present_vectors = str_split(present, ","),
      absent = map(present_vectors, ~ setdiff(accessions_all, .x)),
      absent = map_chr(absent, ~ paste(.x, collapse = ","))
    ) %>%
    select(!present_vectors) %>%
    
    # Calculate MAF: strategy count commas, mutate twice to consider absents with NAs 
    
    ## Notice, absent  = map_chr(absent, ~ paste(.x, collapse = ","))  will create a vector with an empty string c("") if absent is empty
    ## That's why I add the following logical test absent_vect == ""
    
    mutate(present_vect       = str_split(present, pattern = ","),
           present_count      = lengths(present_vect),
           present_count      = if_else(is.na(present_vect), present_count - 1, present_count),
           absent_vect        = str_split(absent, pattern = ","),
           absent_count       = lengths(absent_vect),
           absent_count       = if_else(is.na(absent_vect) | absent_vect == "", absent_count - 1, absent_count),
           total              = present_count + absent_count,
           MAF                = if_else(present_count < absent_count, 
                                        (present_count / total ) * 100, 
                                        (absent_count / total ) * 100 )
    )
  
  final_df <- if (isTRUE(insertions)) {
    intermediate_df %>%
      select(-present_vect, -absent_vect) %>%
      separate(col = TE_id, sep = '-', into = c('chr_rm2', 'start_rm2', 'end_rm2', 'chr_rm', 'start_rm', 'end_rm', 'TE_name'), 
               remove = FALSE ) %>%
      separate(TE_name, sep = ',', into = c('TE_name'), remove = FALSE ) %>%
      separate(TE_name, sep = '#', into = c('TE_numb', 'TE_class'), remove = FALSE ) %>%
      select(!c('chr_rm2', 'start_rm2', 'end_rm2', 'chr_rm', 'start_rm', 'end_rm'))
  } else {
    intermediate_df %>%
      select(-present_vect, -absent_vect) %>%
      separate(col = TE_id, sep = '-', into = c('chr_rm', 'start_rm', 'end_rm', 'TE_name'), remove = FALSE ) %>%
      separate(TE_name, sep = '#', into = c('TE_numb', 'TE_class'), remove = FALSE ) %>%
      select(-chr_rm, -start_rm, -end_rm)
  }
  
  final_df
}