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
        actual <- fix_colors(df)
        expected <- c("#333333", "blue", "#123FFF", "#333333")

        expect_equal(actual, expected)
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
    split_coverages <- data.frame(
        nt_position = c(
            0:4,
            0:9,
            0:3
        ),
        split_name = c(
            rep("contig_1", 5),
            rep("contig_2", 10),
            rep("contig_3", 4)
        ),
        sample_name = rep("sample_1", 19),
        x_values = 0:18
    )

    contig_offsets <- list(
        contig_1 = 0,
        contig_2 = 5, # length of contig 1
        contig_3 = 15 # length of contig 1 + length of contig 2
    )

    describe("get_contig_offsets()", {
        ## This expects data to be properly sorted.
        it("returns list of offsets for each contig", {
            actual <- get_contig_offsets(split_coverages)

            expect_equal(actual, contig_offsets)
        })
    })

    describe("get_offset_snv_positions()", {
        ## Not all contigs need to be in the SNV file.  In the case of
        ## multiple contigs, we need to update the positions in the SNV file
        ## to reflect the "total sample" x-values.

        ## SNV data positions are also 0-based? (TODO CHECK THIS)

        c1_pos <- c(0, 2, 4)
        c2_pos <- c(0, 4, 9)
        c3_pos <- c(0, 1, 3)
        snv_data <- data.frame(
            ## Each contig has 3 SNVs, at the start, end, and somewhere in the
            ## middle.
            contig_name = c(
                rep("contig_1", 3),
                rep("contig_2", 3),
                rep("contig_3", 3)
            ),
            pos_in_contig = c(c1_pos, c2_pos, c3_pos)
        )

        it("returns a thing of offset positions", {
            actual <- get_offset_snv_positions(snv_data, contig_offsets)

            expected <- c(c1_pos, c2_pos + 5, c3_pos + 15)

            expect_equal(actual, expected)
        })
    })



    describe("sort_snv_data()", {

        good_snv_data <- data.frame(
            departure_from_reference = c(0, 1, 0.234)
        )
    })

})
