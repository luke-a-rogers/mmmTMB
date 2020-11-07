# Create mT

# Times: 30
# Areas: 03
# Groups: 01

mT_t30_a03_g01 <- create_release_matrix()

# Write to data/
usethis::use_data(mT_t30_a03_g01, overwrite = TRUE)
