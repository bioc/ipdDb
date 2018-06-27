context("dbMethods")

hla <- loadHlaData()
kir <- loadKirData()
test_that("HLA db can be loaded", {
  expect_is(hla, "IpdDb")
})

test_that("HLA db has valid version", {
  version <- hla$getDbVersion()
  expect_is(version, "character")
  expect_true(grepl("\\d+\\.\\d+\\.\\d+", version))
})

test_that("KIR db can be loaded", {
  expect_is(kir, "IpdDb")
})

test_that("KIR db has valid version", {
  version <- kir$getDbVersion()
  expect_is(version, "character")
  expect_true(grepl("\\d+\\.\\d+\\.\\d+", version))
})

test_that("keytypes and keys method works", {
  kt <- keytypes(hla)
  expect_is(kt, "character")
  ## test that all keytypes are valid and return their keys
  expect_true(all(vapply(kt, function(k, hla) {
    is.character(keys(hla, k))
    }, FUN.VALUE = logical(1), hla = hla)))
  ## Test that a wrong keytype does not work
  expect_error(keys(hla, "testWrongKey"))
})

test_that("columns method works", {
  kt <- keytypes(hla)[1]
  cols <- columns(hla)
  keys <- keys(hla, keytype = kt)[1]
  expect_is(cols, "character")
  ## test that all columns are valid
  expect_true(all(tolower(cols) %in% colnames(select(hla, keys, cols, kt))))
  ## Test if a wrong column gives an error
  expect_error(select(hla, keys, c("blubb", cols), kt))
})

test_that("select method work", {
  kt <- keytypes(hla)[1]
  cols <- columns(hla)[1:2]
  keys <- keys(hla, kt)[1:5]
  res <- select(hla, keys, cols, kt)
  expect_is(res, "data.frame")
  expect_equal(length(res), length(cols))
  expect_equal(length(unique(res[[tolower(kt)]])), length(keys))
})
