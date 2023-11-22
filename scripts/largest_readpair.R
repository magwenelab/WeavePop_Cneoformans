
library(dplyr)
library(readr)

read_csv("read_pair_table.csv") |>
    group_by(sample) |>
    filter(size == max(size, na.rm = TRUE)) |>
    write_csv("largest_read_pair_table.csv")
