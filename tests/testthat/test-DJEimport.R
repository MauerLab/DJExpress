test_that("DJEimport output is class list", {
  in.file <- system.file("extdata", "junct.quant", package = "DJExpress")
  out.file <- DJEimport(in.file, aligner="STAR")
  expect_type(out.file, "list")
})

test_that("DJEimport input contains aligner type", {
  in.file <- system.file("extdata", "junct.quant", package = "DJExpress")
  expect_error(DJEimport(in.file))
})

test_that("DJEimport input contains gtf file path", {
  expect_error(DJEimport(gtf))
})
