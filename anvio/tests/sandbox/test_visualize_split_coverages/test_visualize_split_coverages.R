## Add logger to the global environment.  Many functions assume that
## exists("logger") is TRUE.
assign("logger", make_logger(1), envir = globalenv())

## Stub out the abort() function used by functions that can terminate
## the script.
assign(
  "abort",
  function(...) {
    stop("abort() would've been called")
  },
  envir = globalenv()
)


## Make an alias for describe called context for increased spec
## clarity.
context <- describe

describe("Sample data file", {
  describe("check_samples_against_split_coverages()", {
    ## The point of this function is to remove any samples listed in
    ## the sample data file that aren't also in the split coverages
    ## file.
    split_coverages <- data.frame(
      sample_name = c(
        rep("sample_1", 10),
        rep("sample_2", 10)
      )
    )

    context("some samples are present", {
      sample_data <- data.frame(
        sample_name = paste0("sample_", 1:3),
        sample_color = rep("blue", 3),
        sample_group = c(1, 2, 1)
      )

      it("returns a copy of sample_data without the missing samples", {
        expected <- sample_data[1:2, ]

        actual <- check_samples_against_split_coverages(sample_data, split_coverages)

        expect_equal(actual, expected)
      })
    })

    context("all samples are present", {
      sample_data <- data.frame(
        sample_name = paste0("sample_", 1:2),
        sample_color = rep("blue", 2),
        sample_group = 1:2
      )

      it("returns a copy of sample_data without the missing samples", {
        expected <- sample_data

        actual <- check_samples_against_split_coverages(sample_data, split_coverages)

        expect_equal(actual, expected)
      })
    })

    context("no samples are present", {
      sample_data <- data.frame(
        sample_name = paste0("SAMPLE_", 1:2),
        sample_color = rep("blue", 2),
        sample_group = 1:2
      )

      it("aborts the program", {
        expect_error(check_samples_against_split_coverages(sample_data, split_coverages))
      })
    })
  })

  describe("fix_colors()", {
    df <- data.frame(
      sample_color = c(
        "arst",
        "blue",
        "#123FFF",
        "123FFF"
      )
    )

    it("returns vec with invalid colors as #333333", {
      actual <- fix_colors(df, list(coverage_plot_color = "#333333"))
      expected <- c("#333333", "blue", "#123FFF", "#333333")

      expect_equal(actual, expected)
    })

    context("when sample_color is not specified", {
      sample_data <- data.frame(
        sample_name = paste0("sample_", 1:2),
        group = letters[1:2]
      )

      opts <- data.frame(
        coverage_plot_color = "#333333"
      )

      it("gives a default color for all samples", {
        expected <- rep(opts$coverage_plot_color, times = 2)

        actual <- fix_colors(sample_data, opts)

        expect_equal(actual, expected)
      })
    })
  })

  describe("check_sample_data_headers()", {
    context("with just color info", {
      sample_data <- data.frame(
        sample_name = paste0("sample_", 1:2),
        sample_color = c("blue", "green")
      )

      it("returns NULL and doesn't raise error", {
        expect_null(check_sample_data_headers(sample_data))
      })
    })

    context("with just group info", {
      sample_data <- data.frame(
        sample_name = paste0("sample_", 1:2),
        sample_group = c("a", "b")
      )

      it("returns NULL and doesn't raise error", {
        expect_null(check_sample_data_headers(sample_data))
      })
    })

    context("with color and group info", {
      sample_data <- data.frame(
        sample_name = paste0("sample_", 1:2),
        sample_color = c("blue", "green"),
        sample_group = c("a", "b")
      )

      it("returns NULL and doesn't raise error", {
        expect_null(check_sample_data_headers(sample_data))
      })
    })

    context("with other columns", {
      sample_data <- data.frame(
        sample_name = paste0("sample_", 1:2),
        sample_color = c("blue", "green"),
        sample_group = c("a", "b"),
        apple = 1:2
      )

      it("aborts the program", {
        expect_error(check_sample_data_headers(sample_data))
      })
    })
  })

  describe("check_sample_groups()", {
    default_group_name <- "__anvi__default__sample__group"

    context("when using reserved default group name", {
      sample_data <- data.frame(
        sample_name = "sample_1",
        sample_group = default_group_name
      )

      it("aborts the program", {
        expect_error(check_sample_groups(sample_data))
      })
    })

    context("when no sample_group info is given", {
      sample_data <- data.frame(
        sample_name = "sample_1",
        sample_color = "#333333"
      )

      it("adds the default group name", {
        expected <- default_group_name
        actual <- check_sample_groups(sample_data)

        expect_equal(actual, expected)
      })
    })
  })
})

describe("sapply_by_group()", {
  describe("with character groups", {
    df <- data.frame(
      count = c(1, 20, 3, 10, 2, 30),
      group = rep(c("a", "b"), 3)
    )

    it("sapply a function by group", {
      actual <- sapply_by_group(count ~ group, df, mean)
      expect_equal(actual, c(2, 20))
    })

    it("follows sort order of the factors", {
      ## Change sort order
      df$group <- factor(df$group, levels = c("b", "a"))

      actual <- sapply_by_group(count ~ group, df, mean)
      expect_equal(actual, c(20, 2))
    })
  })

  describe("with numeric groups", {
    df <- data.frame(
      count = c(1, 20, 3, 10, 2, 30),
      group = rep(c(1, 2), 3)
    )

    it("sapply a function by group", {
      actual <- sapply_by_group(count ~ group, df, mean)
      expect_equal(actual, c(2, 20))
    })
  })
})

describe("get_points_per_sample()", {
  total <- 10

  df <- data.frame(
    coverage = 1:10,
    split_name = rep("contig", total),
    sample_name = c(rep("sample_1", total / 2),
                    rep("sample_2", total / 2))
  )

  it("calculates number of points per sample", {
    actual <- get_points_per_sample(df)
    expected <- rep(total / 2, 2)

    expect_equal(actual, expected)
  })
})

describe("smooth_coverage()", {
  coverage <- c(1, 4, 3, 4)
  smoothed_coverage <- c(1, 3, 4, 4)

  df <- data.frame(
    coverage = c(coverage,
                 coverage * 10,
                 coverage * 100,
                 coverage * 1000),
    split_name = c(rep("contig_1", 4),
                   rep("contig_2", 4),
                   rep("contig_1", 4),
                   rep("contig_2", 4)),
    ## sample_2 comes before sample_1 in the data frame!
    sample_name = c(rep("sample_2", 8),
                    rep("sample_1", 8))
  )

  it("smooths coverage with running medians", {
    actual <- smooth_coverage(df, window_size = 3)

    ## sample_1 values come first in the output however as they
    ## are grouped by sample, and that is sorted by factor level!
    expected <- c(
      smoothed_coverage * 100,
      smoothed_coverage * 1000,
      smoothed_coverage,
      smoothed_coverage * 10
    )

    expect_equal(actual, expected)
  })
})

describe("shrink_data()", {
  df <- data.frame(
    sample_name = c(rep("s1", 8),
                    rep("s2", 8)),
    x_values = rep(0:7, 2),
    split_name = c(rep("c1", 4),
                   rep("c2", 4),
                   rep("c1", 4),
                   rep("c2", 4)),
    coverage = rep(10, 16)
  )

  it("reduces the size of the data frame", {
    ## For each sample, keeps first and last rows plus any that
    ## divide the window size.
    actual <- shrink_data(df, window_size = 3)
    expected <- df[c(1, 4, 7, 8, 9, 12, 15, 16), ]
    expected$x_values.max <- rep(7, 8)

    expect_equal(actual, expected)
  })
})

describe("is_hex_color()", {
  it("hex codes must start with '#' char", {
      expect_false(is_hex_color("123fff"))
  })

  it("is false for 'short' hex codes", {
    expect_false(is_hex_color("#333"))
  })

  it("is true for valid hex code", {
    expect_true(is_hex_color("#123fFf"))
    expect_true(is_hex_color("#000000"))
    expect_true(is_hex_color("#fFfFfF"))
  })

  it("is false for invalid hex code", {
    expect_false(is_hex_color("#123qqq"))
    expect_false(is_hex_color("black"))
    expect_false(is_hex_color("#123q"))
  })
})


describe("sort_split_coverages()", {
  sample_sort_order <- c("s_2", "s_1")

  split_coverages <- data.frame(
    unique_entry_id = 0:9,
    nt_position = rep(0:4, 2),
    split_name = rep("contig_1", 10),
    sample_name = c(rep("s_1", 5), rep("s_2", 5)),
    coverage = 1:10,
    x_values = rep(0:4, 2)
  )

  sample_data <- data.frame(
    sample_name = sample_sort_order,
    sample_color = c("blue", "green")
  )

  it("sorts split_coverages df based on sample ordering in sample_data", {
    actual <- sort_split_coverages(split_coverages, sample_data)
    expected <- rbind(
      split_coverages[6:10, ],
      split_coverages[1:5, ]
    )

    expect_equivalent(actual, expected)
  })
})

## SNV data functions
describe("SNV data function", {
  expected_contig_offsets <- matrix(
    rep(c(0, 5, 15), 2),
    nrow = 3,
    ncol = 2,
    dimnames = list(
      split_name = c("contig_1_split_00001",
                     "contig_1_split_00002",
                     "contig_2_split_00001"),
      sample_name = paste0("sample_", 1:2)
    )
  )


  describe("get_contig_offsets()", {
    describe("with multiple samples and multiple contigs", {
      split_coverages <- data.frame(
        nt_position = rep(
          c(0:4, 0:9, 0:3),
          times = 2
        ),
        split_name = rep(
          c(rep("contig_1_split_00001", 5),
            rep("contig_1_split_00002", 10),
            rep("contig_2_split_00001", 4)),
          times = 2
        ),
        sample_name = c(
          rep("sample_1", 19),
          rep("sample_2", 19)
        ),
        x_values = rep(
          0:18,
          times = 2
        )
      )

      it("returns a data.frame with length of contigs in samples", {
        actual_contig_offsets <- get_contig_offsets(split_coverages)

        expect_equal(actual_contig_offsets, expected_contig_offsets)
      })
    })

    describe("with multiple samples and a single contig", {
      split_coverages <- data.frame(
        nt_position = rep(0:4, 2),
        split_name = rep("contig_1_split_00001", 10),
        sample_name = c(rep("sample_1", 5), rep("sample_2", 5)),
        x_values = rep(0:4, 2)
      )

      expected_contig_offsets <- matrix(
        c(0, 0),
        nrow = 1,
        ncol = 2,
        dimnames = list(
          split_name = c("contig_1_split_00001"),
          sample_name = paste0("sample_", 1:2)
        )
      )

      it("returns a 1 x nsamples data frame", {
        actual_contig_offsets <- get_contig_offsets(split_coverages)

        expect_equal(actual_contig_offsets, expected_contig_offsets)
      })
    })

    describe("with a single sample and multiple contigs", {
      split_coverages <- data.frame(
        nt_position = c(0:4, 0:9),
        split_name = c(
          rep("contig_1_split_00001", 5),
          rep("contig_1_split_00002", 10)
        ),
        sample_name = rep("sample_1", 15),
        x_values = 0:14
      )

      expected_contig_offsets <- matrix(
        c(0, 5),
        nrow = 2,
        ncol = 1,
        dimnames = list(
          split_name = c("contig_1_split_00001", "contig_1_split_00002"),
          sample_name = c("sample_1")
        )
      )

      it("returns a 1 x nsamples data frame", {
        actual_contig_offsets <- get_contig_offsets(split_coverages)

        expect_equal(actual_contig_offsets, expected_contig_offsets)
      })
    })


    describe("with single sample and a single contig", {
      split_coverages <- data.frame(
        nt_position = 0:4,
        split_name = rep("contig_1_split_00001", 5),
        sample_name = rep("sample_1", 5),
        x_values = 0:4
      )

      expected_contig_offsets <- matrix(
        0,
        nrow = 1,
        ncol = 1,
        dimnames = list(
          split_name = c("contig_1_split_00001"),
          sample_name = c("sample_1")
        )
      )

      it("returns a 1 x nsamples data frame", {
        actual_contig_offsets <- get_contig_offsets(split_coverages)

        expect_equal(actual_contig_offsets, expected_contig_offsets)
      })
    })

  })

  describe("get_offset_snv_positions()", {
    ## Not all contigs need to be in the SNV file.  In the case of
    ## multiple contigs, we need to update the positions in the SNV file
    ## to reflect the "total sample" x-values.

    ## Technically, this should be sorted on sample, then by contig,
    ## but it doesn't matter for the test.
    snv_data <- data.frame(
      ## Each contig has 3 SNVs, at the start, end, and somewhere in the
      ## middle.
      split_name = c(
        ## each has 3 SNVs and in 2 samples
        rep("contig_1_split_00001", 3 * 2),
        rep("contig_1_split_00002", 3 * 2),
        rep("contig_2_split_00001", 3 * 2)
      ),
      sample_name = rep(paste0("sample_", 1:2), 3),
      pos = c(
        ## Contig 1
        0, 0,
        2, 2,
        4, 4,
        ## Contig 2
        0, 0,
        4, 4,
        9, 9,
        ## Contig 3 snv positions
        0, 0,
        1, 1,
        3, 3
      )
    )

    it("returns a vector of offset positions", {
      actual <- get_offset_snv_positions(snv_data, expected_contig_offsets)

      c1_offset <- 0
      c2_offset <- 5
      c3_offset <- 15

      expected <- c(
        ## Contig 1
        0 + c1_offset, 0 + c1_offset,
        2 + c1_offset, 2 + c1_offset,
        4 + c1_offset, 4 + c1_offset,

        ## Contig 2
        0 + c2_offset, 0 + c2_offset,
        4 + c2_offset, 4 + c2_offset,
        9 + c2_offset, 9 + c2_offset,

        ## Contig 3
        0 + c3_offset, 0 + c3_offset,
        1 + c3_offset, 1 + c3_offset,
        3 + c3_offset, 3 + c3_offset
      )

      expect_equal(actual, expected)
    })
  })
})
