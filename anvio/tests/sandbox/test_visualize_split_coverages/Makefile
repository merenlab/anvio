TEST_FILE_DIR = test_files
TEST_OUT_DIR = TEST_OUTDIR
PROJ_ROOT_DIR = ../../../..

SCRIPT = $(PROJ_ROOT_DIR)/sandbox/anvi-script-visualize-split-coverages

.DEFAULT_GOAL := test_all

test_full:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output && ls -lah $(TEST_OUT_DIR) && echo

test_full_line:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output --chart-type line && ls -lah $(TEST_OUT_DIR) && echo

test_no_color:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data_no_color.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output && ls -lah $(TEST_OUT_DIR)  && echo

# There was a bug with xlim, area, and no sample colors.
test_no_color_area_xlim:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data_no_color.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output --chart-type area --xlim 10,20 && ls -lah $(TEST_OUT_DIR)  && echo

test_no_sample_group:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data_no_group.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output && ls -lah $(TEST_OUT_DIR)  && echo

test_no_sample_data:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output && ls -lah $(TEST_OUT_DIR) && echo

# There was a bug with xlim, area, and no sample colors.
test_no_sample_data_area_xlim:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output --chart-type area --xlim 10,20 && ls -lah $(TEST_OUT_DIR) && echo

test_no_snv_data:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output && ls -lah $(TEST_OUT_DIR) && echo

test_no_sample_data_no_snv_data:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output && ls -lah $(TEST_OUT_DIR) && echo

test_with_pdf_extension:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output.pdf && ls -lah $(TEST_OUT_DIR) && echo

test_with_pdf_extension_no_sample_group:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data_no_group.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output.pdf && ls -lah $(TEST_OUT_DIR) && echo

test_xlim:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output --xlim 10,20 && ls -lah $(TEST_OUT_DIR) && echo

test_xlim_area:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output --xlim 10,20 --chart-type area && ls -lah $(TEST_OUT_DIR) && echo

test_xlim_line:
	rm -r $(TEST_OUT_DIR); mkdir -p $(TEST_OUT_DIR); $(SCRIPT) --infile $(TEST_FILE_DIR)/split_cov.txt --sample-data $(TEST_FILE_DIR)/sample_data.txt --snv-data $(TEST_FILE_DIR)/snv.txt --outfile-basename $(TEST_OUT_DIR)/vis_split_cov_output --xlim 10,20 --chart-type line && ls -lah $(TEST_OUT_DIR) && echo

test_all: test_full test_full_line test_no_color test_no_color_area_xlim test_no_sample_group test_no_sample_data test_no_sample_data_area_xlim test_no_snv_data test_no_sample_data_no_snv_data test_with_pdf_extension test_with_pdf_extension_no_sample_group test_xlim test_xlim_area test_xlim_line
