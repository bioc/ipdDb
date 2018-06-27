context("seqMethods")
hla <- loadHlaData()

test_that("getReference works as expected", {
  locus <- getLoci(hla)[1]
  alleles <- getAlleles(hla, locus)
  ## Test single allele
  seq <- getReference(hla, alleles[1])
  expect_is(seq, "DNAStringSet")
  ## Test multiple alleles
  seq <- getReference(hla, alleles[1:10])
  expect_is(seq, "DNAStringSet")
  expect_error(getReference(hla, "xx"))
  expect_error(getReference(hla, 1))
  expect_error(getReference(hla, TRUE))
  expect_error(getReference(hla, locus))
})

test_that("getLoci works as expected", {
  loci <- getLoci(hla)
  expect_is(loci, "character")
  ## check if all loci begin with "HLA" or "KIR"
  expect_true(all(vapply(loci, function(l) startsWith(l, "HLA-") || 
                           startsWith(l, "KIR"), FUN.VALUE = logical(1))))
})

test_that(".matchAllele finds correct alleles", {
  locus <- "HLA-A"
  allAlleles <- getAlleles(hla, locus)
  targetAllele <- "HLA-A*01:01:01:01"
  
  ## test partial alleles 
  allele <- "01:01:01:01"
  expect_identical(.matchAllele(hla, allele, locus), targetAllele)
  
  ## test truncated alleles; shorter resolution
  allele <- "HLA-A*01:01:01"
  targetAlleles <- allAlleles[startsWith(allAlleles, allele)]
  expect_identical(.matchAllele(hla, allele, locus), targetAlleles)
  
  allele <- "HLA-A*9999"
  expect_identical(.matchAllele(hla, allele, locus), allAlleles)
  
  expect_error(.matchAllele(hla, 123, locus))
  expect_error(.matchAllele(hla, TRUE, locus))
  expect_error(.matchAllele(hla, allele, "testlocus"))
})

test_that("getAlleles works correct", {
  locus <- "HLA-A"
  ## Get the target alleles using plain methods
  keys <- keys(hla, keytype = "locus" )
  targetAlleles <- select(hla, keys = locus, 
                          columns = "allele_name", 
                          keytype = "locus")$allele_name
  expect_identical(getAlleles(hla, locus), targetAlleles)
  expect_error(getAlleles(hla, "invalidLocus"))
})

test_that("getStructure works correct",{
  alleles <- getAlleles(hla, getLoci(hla)[1])
  expect_is(getStructure(hla, alleles[1]), "GRanges")
  expect_is(getStructure(hla, alleles[1:10]), "GRanges")
  
  ## Test for invalid keys
  expect_error(getStructure(hla, "test"))
  expect_error(getStructure(hla, 1))
  expect_error(getStructure(hla, TRUE))
})

test_that("closestComplete gives correct sequences", {
  alleles <- getAlleles(hla, getLoci(hla)[1])
  expect_is(getClosestComplete(hla, alleles[1]), "DNAStringSet")
  expect_error(getClosestComplete(hla, alleles[1:100]))
  expect_error(getClosestComplete(hla, "test"))
  expect_error(getClosestComplete(hla, TRUE))
  expect_error(getClosestComplete(hla, 123))
})
