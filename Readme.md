ipdDb
=============

Database package of HLA and KIR alleles from the 
[IPD IMGT/HLA](https://www.ebi.ac.uk/ipd/imgt/hla/) and 
[IPD KIR](https://www.ebi.ac.uk/ipd/kir/) database.
Reference:

>Robinson J, Maccari G, Marsh SGE, Walter L, Blokhuis J, Bimber B, Parham P, De Groot NG, Bontrop RE, Guethlein LA, and Hammond JA
>KIR Nomenclature in non-human species
>Immunogenetics (2018), in preparation



This package holds all information from the IPD databases.
For HLA this is limited to HLA-A, -B, -C, -DPB1, -DQB1 and -DRB1.

For alleles which are not known in full-length, also the closest full-length allele is stored.

The data is stored in an SQLite database and its compatible with Bioconductor's database package rules, 
i.e. the `select` method including `columns`, `keytypes` and `keys`.

Additionally, some helper functions are implemented to fetch all contained loci, all alleles of a locus,
the sequences of alleles and the sequence of closest allele which is available in full-length for an allele.

```r
library(ipdDb)

hla <- loadHlaData()
## get all loci stored in the db
available_loci <- hla$getLoci()

## get all alleles of a locus
alleles <- hla$getAlleles(available_loci[1]))
alleles <- hla$getAlleles("HLA-A")

## get all sequences of a bunch of alleles as DNAStringSet
sequences <- hla$getReferences(alleles)
sequences <- hla$getReferences(c("HLA-A*01:01:01:01", "HLA-A*01:01:01:03" ))

## get the closest complete reference for ONE allele as DNAStringSet
closest_complete <- hla$getClosestComplete(alleles[1])
closest_complete <- hla$getClosestComplete("HLA-A*01:01:01:01")

## Get the gene structure for a bunch of alleles as GRanges object
structures <- hla$getStructure(ipd.Hsapiens.db, alleles)
structures <- hla$getStructure(ipd.Hsapiens.db, c("HLA-A*01:01:01:01", "HLA-A*01:01:01:03" ))

```


