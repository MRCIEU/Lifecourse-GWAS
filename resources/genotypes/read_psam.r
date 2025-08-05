library(data.table)

#' Read and validate PSAM (PLINK Sample) files
#'
#' @param file Path to the PSAM file
#' @param validate Logical, whether to perform validation (default: TRUE)
#' @return A data.table with the PSAM data and attributes for metadata
#' @export
read_psam <- function(file, validate = TRUE) {
  
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }
  
  # Read all lines to parse header and body
  lines <- readLines(file, warn = FALSE)
  
  if (length(lines) == 0) {
    stop("File is empty")
  }
  
  # Identify header lines (start with #)
  header_lines <- grep("^#", lines)
  
  # Parse header
  if (length(header_lines) == 0) {
    # No header present - use default based on number of columns
    first_data_line <- lines[1]
    # Split on tabs or spaces
    fields <- strsplit(first_data_line, "[ \t]+")[[1]]
    n_cols <- length(fields)
    
    if (n_cols >= 6) {
      header <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENO1")
    } else if (n_cols == 5) {
      header <- c("FID", "IID", "PAT", "MAT", "SEX")
    } else {
      stop("Invalid number of columns for headerless file: ", n_cols)
    }
    
    body_start <- 1
  } else {
    # Find the last header line that starts with #FID or #IID
    column_header_idx <- NULL
    for (i in rev(header_lines)) {
      line <- lines[i]
      if (grepl("^#(FID|IID)", line)) {
        column_header_idx <- i
        break
      }
    }
    
    if (is.null(column_header_idx)) {
      stop("No valid column header line found (must start with #FID or #IID)")
    }
    
    # Parse column headers
    header_line <- lines[column_header_idx]
    header <- strsplit(gsub("^#", "", header_line), "[ \t]+")[[1]]
    
    # Remove empty strings
    header <- header[header != ""]
    
    body_start <- max(header_lines) + 1
  }
  
  # Validate header structure
  if (validate) {
    validate_header(header)
  }
  
  # Read body data
  if (body_start > length(lines)) {
    stop("No data rows found in file")
  }
  
  body_lines <- lines[body_start:length(lines)]
  # Remove empty lines
  body_lines <- body_lines[body_lines != ""]
  
  if (length(body_lines) == 0) {
    stop("No non-empty data rows found")
  }
  
  # Parse body data
  data_list <- lapply(body_lines, function(line) {
    strsplit(line, "[ \t]+")[[1]]
  })
  
  # Check that all rows have at least as many columns as header
  min_cols <- length(header)
  row_lengths <- sapply(data_list, length)
  
  if (any(row_lengths < min_cols)) {
    bad_rows <- which(row_lengths < min_cols)
    stop("Rows with insufficient columns found at lines: ", 
         paste(body_start + bad_rows - 1, collapse = ", "))
  }
  
  # Truncate rows to header length (ignore extra columns)
  data_list <- lapply(data_list, function(row) row[1:min_cols])
  
  # Convert to data.table
  dt <- data.table(do.call(rbind, data_list))
  setnames(dt, header)
  
  # Perform validation if requested
  if (validate) {
    validate_psam_data(dt)
  }
  
  # Add metadata as attributes
  attr(dt, "header") <- header
  attr(dt, "n_header_lines") <- length(header_lines)
  
  return(dt)
}

#' Validate PSAM header structure
#'
#' @param header Character vector of column names
validate_header <- function(header) {
  
  # Check for duplicate headers
  if (any(duplicated(header))) {
    stop("Duplicate column headers not allowed")
  }
  
  # Check FID/IID positioning
  fid_pos <- match("FID", header)
  iid_pos <- match("IID", header)
  
  if (is.na(iid_pos)) {
    stop("IID column is required")
  }
  
  if (!is.na(fid_pos)) {
    if (fid_pos != 1) {
      stop("FID column must be first if present")
    }
    if (iid_pos != 2) {
      stop("IID column must immediately follow FID column")
    }
  } else {
    if (iid_pos != 1) {
      stop("IID column must be first if FID column is absent")
    }
  }
  
  # Check SID positioning
  sid_pos <- match("SID", header)
  if (!is.na(sid_pos)) {
    expected_sid_pos <- if (is.na(fid_pos)) 2 else 3
    if (sid_pos != expected_sid_pos) {
      stop("SID column must immediately follow IID column")
    }
  }
  
  # Check PAT/MAT presence
  pat_pos <- match("PAT", header)
  mat_pos <- match("MAT", header)
  
  if (is.na(pat_pos) != is.na(mat_pos)) {
    stop("PAT and MAT columns must both be present or both be absent")
  }
}

#' Validate PSAM data content
#'
#' @param dt data.table with PSAM data
validate_psam_data <- function(dt) {
  
  # Check IID values are not '0'
  if ("IID" %in% names(dt)) {
    if (any(dt$IID == "0")) {
      stop("IID values cannot be '0'")
    }
  }
  
  # Check sample ID uniqueness
  sample_ids <- create_sample_ids(dt)
  if (any(duplicated(sample_ids))) {
    dup_ids <- sample_ids[duplicated(sample_ids)]
    stop("Duplicate sample IDs found: ", paste(unique(dup_ids), collapse = ", "))
  }
  
  # Validate SEX column if present
  if ("SEX" %in% names(dt)) {
    valid_sex <- c("1", "2", "M", "m", "F", "f", "NA", "0", "-9")
    invalid_sex <- dt$SEX[!dt$SEX %in% valid_sex]
    if (length(invalid_sex) > 0) {
      warning("Invalid SEX values found (will be treated as missing): ", 
              paste(unique(invalid_sex), collapse = ", "))
    }
  }
  
  # Check phenotype columns for proper formatting
  predefined_cols <- c("FID", "IID", "SID", "PAT", "MAT", "SEX")
  pheno_cols <- setdiff(names(dt), predefined_cols)
  
  for (col in pheno_cols) {
    validate_phenotype_column(dt[[col]], col)
  }
}

#' Create full sample IDs (FID-IID-SID)
#'
#' @param dt data.table with PSAM data
#' @return Character vector of sample IDs
create_sample_ids <- function(dt) {
  fid <- if ("FID" %in% names(dt)) dt$FID else rep("0", nrow(dt))
  iid <- dt$IID
  sid <- if ("SID" %in% names(dt)) dt$SID else rep("0", nrow(dt))
  
  paste(fid, iid, sid, sep = "-")
}

#' Validate phenotype column
#'
#' @param values Character vector of phenotype values
#' @param col_name Name of the column for error reporting
validate_phenotype_column <- function(values, col_name) {
  
  # Check if any value starts with digit (indicating quantitative/binary)
  starts_with_digit <- grepl("^[+-]?[0-9]", values)
  
  # Check for NA variants
  na_variants <- toupper(values) %in% c("NA", "NAN")
  
  if (any(starts_with_digit) || any(na_variants)) {
    # This appears to be a quantitative or binary phenotype
    # Check for binary phenotype pattern
    non_missing <- values[!na_variants & values != "-9" & values != "0"]
    
    if (all(non_missing %in% c("1", "2"))) {
      # This looks like a binary phenotype - that's fine
      return(invisible(TRUE))
    }
    
    # Check if all values are numeric
    numeric_values <- suppressWarnings(as.numeric(values[!na_variants]))
    if (any(is.na(numeric_values))) {
      warning("Phenotype column '", col_name, 
              "' contains non-numeric values but starts with digits")
    }
  } else {
    # This appears to be a categorical phenotype
    # Check that no values start with digits
    if (any(starts_with_digit)) {
      stop("Categorical phenotype column '", col_name, 
           "' contains values starting with digits")
    }
    
    # Check for invalid missing value representations
    invalid_missing <- toupper(values) %in% c("NA", "NAN")
    if (any(invalid_missing)) {
      warning("Categorical phenotype column '", col_name, 
              "' uses 'NA'/'NAN' for missing values. Use 'NONE' instead.")
    }
  }
  
  invisible(TRUE)
}

# Unit Tests
run_psam_tests <- function() {
  
  cat("Running PSAM reader unit tests...\n\n")
  
  # Test 1: Basic PSAM file with header
  cat("Test 1: Basic PSAM file with header\n")
  test_content1 <- c(
    "#FID IID PAT MAT SEX PHENO1",
    "FAM1 ID1 0 0 1 1",
    "FAM1 ID2 0 0 2 2",
    "FAM2 ID3 ID1 ID2 1 1"
  )
  
  temp_file1 <- tempfile(fileext = ".psam")
  writeLines(test_content1, temp_file1)
  
  tryCatch({
    dt1 <- read_psam(temp_file1)
    print(str(dt1))
    stopifnot(nrow(dt1) == 3)
    stopifnot(ncol(dt1) == 6)
    stopifnot(all(names(dt1) == c("FID", "IID", "PAT", "MAT", "SEX", "PHENO1")))
    cat("✓ PASSED\n\n")
  }, error = function(e) {
    cat("✗ FAILED:", e$message, "\n\n")
  })
  
  unlink(temp_file1)
  
  # Test 2: Headerless file with 6 columns
  cat("Test 2: Headerless file with 6 columns\n")
  test_content2 <- c(
    "FAM1 ID1 0 0 1 1.5",
    "FAM1 ID2 0 0 2 2.3"
  )
  
  temp_file2 <- tempfile(fileext = ".psam")
  writeLines(test_content2, temp_file2)
  
  tryCatch({
    dt2 <- read_psam(temp_file2)
    print(str(dt2))
    stopifnot(nrow(dt2) == 2)
    stopifnot(all(names(dt2) == c("FID", "IID", "PAT", "MAT", "SEX", "PHENO1")))
    cat("✓ PASSED\n\n")
  }, error = function(e) {
    cat("✗ FAILED:", e$message, "\n\n")
  })
  
  unlink(temp_file2)
  
  # Test 3: File with SID column
  cat("Test 3: File with SID column\n")
  test_content3 <- c(
    "#FID IID SID PAT MAT SEX",
    "FAM1 ID1 S1 0 0 1",
    "FAM1 ID1 S2 0 0 1"
  )
  
  temp_file3 <- tempfile(fileext = ".psam")
  writeLines(test_content3, temp_file3)
  
  tryCatch({
    dt3 <- read_psam(temp_file3)
    print(str(dt3))
    stopifnot(nrow(dt3) == 2)
    stopifnot("SID" %in% names(dt3))
    cat("✓ PASSED\n\n")
  }, error = function(e) {
    cat("✗ FAILED:", e$message, "\n\n")
  })
  
  unlink(temp_file3)
  
  # Test 4: Invalid IID values (should fail)
  cat("Test 4: Invalid IID values (should fail)\n")
  test_content4 <- c(
    "#FID IID PAT MAT SEX",
    "FAM1 0 0 0 1"  # IID cannot be 0
  )
  
  temp_file4 <- tempfile(fileext = ".psam")
  writeLines(test_content4, temp_file4)
  
  tryCatch({
    dt4 <- read_psam(temp_file4)
    print(str(dt4))
    cat("✗ FAILED: Should have thrown error for IID = 0\n\n")
  }, error = function(e) {
    if (grepl("IID values cannot be '0'", e$message)) {
      cat("✓ PASSED (correctly caught invalid IID)\n\n")
    } else {
      cat("✗ FAILED: Wrong error message:", e$message, "\n\n")
    }
  })
  
  unlink(temp_file4)
  
  # Test 5: Duplicate sample IDs (should fail)
  cat("Test 5: Duplicate sample IDs (should fail)\n")
  test_content5 <- c(
    "#FID IID PAT MAT SEX",
    "FAM1 ID1 0 0 1",
    "FAM1 ID1 0 0 2"  # Duplicate sample ID
  )
  
  temp_file5 <- tempfile(fileext = ".psam")
  writeLines(test_content5, temp_file5)
  
  tryCatch({
    dt5 <- read_psam(temp_file5)
    print(str(dt5))
    cat("✗ FAILED: Should have thrown error for duplicate sample IDs\n\n")
  }, error = function(e) {
    if (grepl("Duplicate sample IDs", e$message)) {
      cat("✓ PASSED (correctly caught duplicate sample IDs)\n\n")
    } else {
      cat("✗ FAILED: Wrong error message:", e$message, "\n\n")
    }
  })
  
  unlink(temp_file5)
  
  # Test 6: Missing PAT without MAT (should fail)
  cat("Test 6: PAT without MAT (should fail)\n")
  test_content6 <- c(
    "#FID IID PAT SEX",  # Missing MAT
    "FAM1 ID1 0 1"
  )
  
  temp_file6 <- tempfile(fileext = ".psam")
  writeLines(test_content6, temp_file6)
  
  tryCatch({
    dt6 <- read_psam(temp_file6)
    print(str(dt6))
    cat("✗ FAILED: Should have thrown error for PAT without MAT\n\n")
  }, error = function(e) {
    if (grepl("PAT and MAT columns must both be present", e$message)) {
      cat("✓ PASSED (correctly caught PAT without MAT)\n\n")
    } else {
      cat("✗ FAILED: Wrong error message:", e$message, "\n\n")
    }
  })
  
  unlink(temp_file6)
  
  # Test 7: Categorical phenotype
  cat("Test 7: Categorical phenotype\n")
  test_content7 <- c(
    "#IID CATEGORY",
    "ID1 TYPE_A",
    "ID2 TYPE_B",
    "ID3 NONE"
  )
  
  temp_file7 <- tempfile(fileext = ".psam")
  writeLines(test_content7, temp_file7)
  
  tryCatch({
    dt7 <- read_psam(temp_file7)
    print(str(dt7))
    stopifnot(nrow(dt7) == 3)
    stopifnot(all(names(dt7) == c("IID", "CATEGORY")))
    cat("✓ PASSED\n\n")
  }, error = function(e) {
    cat("✗ FAILED:", e$message, "\n\n")
  })
  
  unlink(temp_file7)

  # Test 8: File with SID column
  cat("Test 8: File with extra headers\n")
  test_content8 <- c(
    "# Extra header line",
    "#FID IID SID PAT MAT SEX",
    "FAM1 ID1 S1 0 0 1",
    "FAM1 ID1 S2 0 0 1"
  )
  
  temp_file8 <- tempfile(fileext = ".psam")
  writeLines(test_content8, temp_file8)
  
  tryCatch({
    dt8 <- read_psam(temp_file8)
    print(str(dt8))
    stopifnot(nrow(dt8) == 2)
    stopifnot("SID" %in% names(dt8))
    cat("✓ PASSED\n\n")
  }, error = function(e) {
    cat("✗ FAILED:", e$message, "\n\n")
  })
  
  unlink(temp_file8)


  
  cat("All tests completed!\n")
}
