library(testthat)

## First get the directory of this test runner file.
opts <- commandArgs()
file_opt_idx <- grep("--file=", opts, fixed = TRUE)
file_opt <- opts[[file_opt_idx]]
file_name <- sub("--file=", "", file_opt, fixed = TRUE)
this_dir <- dirname(file_name)

source_code <- file.path(
    this_dir,
    "..", "..", "..", "..",
    "sandbox",
    "anvi-script-visualize-split-coverages"
)

## If this line is changed in the source file, this hack won't work
## anymore.
function_end_marker <- "End of function definitions"
source_code_lines <- readLines(source_code)
function_end_line <- grep(function_end_marker, source_code_lines, fixed = TRUE)

## Only source the part of the file with function definitions.
tmp <- textConnection(source_code_lines[1:function_end_line])
source(tmp)
close(tmp)

## Run the tests.
test_dir(this_dir, reporter = SummaryReporter$new())
