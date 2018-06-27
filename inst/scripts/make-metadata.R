hlaMetadata <- data.frame(
  Title = "Allele data from the IPD IMGT/HLA database",
  Description = paste0("Data for all alleles of selected HLA loci (HLA-A, -B, 
  -C, -DPB1, -DQB1 and -DRB1). The allele annotation, sequence, gene structure
  and the (sequence-based) closest allele in full-length is stored. 
  Reference:
  Robinson J, Maccari G, Marsh SGE, Walter L, Blokhuis J, Bimber B, Parham P, 
  De Groot NG, Bontrop RE, Guethlein LA, and Hammond JA
  KIR Nomenclature in non-human species Immunogenetics (2018), in preparation"),
  BiocVersion = "3.8",
  Genome = "no genome data",
  SourceType = "Zip",
  SourceUrl = "https://github.com/ANHIG/IMGTHLA/blob/Latest/xml/hla.xml.zip",
  SourceVersion = "3.32.0",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = "EMBL-EBI",
  Maintainer = "Steffen Klasberg <klasberg@dkms-lab.de>",
  RDataClass = "data.frame, DNAStringSet, GRanges",
  DispatchClass = "SQLiteFile",
  RDataPath = "ipdDb/ipdHLA_3.32.0.sqlite",
  Tags = "ipd:hla:IMGT/HLA:alleles"
)
kirMetadata <- data.frame(
  Title = "Allele data from the IPD KIR database",
  Description = paste0("Data for the alleles of all KIR loci in the database. 
  The allele annotation, sequence, gene structure and the (sequence-based) 
  closest allele in full-length is stored. 
  Reference:
  Robinson J, Maccari G, Marsh SGE, Walter L, Blokhuis J, Bimber B, Parham P, 
  De Groot NG, Bontrop RE, Guethlein LA, and Hammond JA
  KIR Nomenclature in non-human species Immunogenetics (2018), in preparation"),
  BiocVersion = "3.8",
  Genome = "no genome data",
  SourceType = "Zip",
  SourceUrl = "https://github.com/ANHIG/IPDKIR/blob/Latest/KIR.dat",
  SourceVersion = "2.7.1",
  Species = "Homo sapiens",
  TaxonomyId = 9606,
  Coordinate_1_based = TRUE,
  DataProvider = "EMBL-EBI",
  Maintainer = "Steffen Klasberg <klasberg@dkms-lab.de>",
  RDataClass = "data.frame, DNAStringSet, GRanges",
  DispatchClass = "SQLiteFile",
  RDataPath = "ipdDb/ipdKIR_2.7.1.sqlite",
  Tags = "ipd:kir:alleles"
)

write.csv(rbind(hlaMetadata, kirMetadata), "inst/extdata/ipd_metadata.csv")
