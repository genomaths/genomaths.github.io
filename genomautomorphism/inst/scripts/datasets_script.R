# ========================================================================= #
#
# ======== Script used to generate the datasets used in the examples ====== #
#
# ========================================================================= #
library(GenomAutomorphism)
library(Biostrings)


aln <- c(
    "ACCTATGTTGGTATT---GCGCTCCAACTCCTTGGCTCTAGCTCACTACAT",
    "ATCTATGTTGGTATTACGACGCTCCAATTCCTTGGGTCC------CTCCTT"
)

aln <- DNAStringSet(aln)

usethis::use_data(aln)

URL <- paste0(
    "https://github.com/genomaths/seqalignments/raw/master/CYCS/",
    "primate_cytochrome_c_(CYCS)_18_sequences.fasta"
)

cyc_aln <- readDNAMultipleAlignment(filepath = URL)

usethis::use_data(cyc_aln, overwrite = TRUE)


nams <- c(
    "human_1", "human_2", "gorilla", "human_3", "human_4",
    "human_5", "human_6", "silvery_gibbon", "white_cheeked_gibbon",
    "francois_langur", "olive_baboon_1", "olive_baboon_2",
    "golden_monkey", "rhesus_monkeys_1", "rhesus_monkeys_2",
    "gelada_baboon_1", "gelada_baboon_2", "orangutan_1", "orangutan_2"
)

cyc_autm <- automorphisms(
    filepath = URL,
    group = "Z64",
    cube = c("ACGT", "TGCA"),
    cube_alt = c("CATG", "GTAC"),
    nms = nams
)


usethis::use_data(cyc_autm, overwrite = TRUE)


URL <- paste0(
    "https://github.com/genomaths/seqalignments/raw/master/BRCA1/",
    "brca1_primates_dna_repair_20_sequences.fasta"
)


nams <- c(
    "human_1", "human_2", "gorilla_1", "gorilla_2", "gorilla_3",
    "chimpanzee_1", "chimpanzee_2", "chimpanzee_3", "chimpanzee_4",
    "bonobos_1", "bonobos_2", "bonobos_3", "bonobos_4", "silvery_gibbon_1",
    "silvery_gibbon_1", "silvery_gibbon_3", "golden_monkey_1",
    "golden_monkey_2", "gelada_baboon", "bolivian_monkey"
)

brca1_aln <- readDNAMultipleAlignment(filepath = URL)


brca1_autm <- automorphisms(
    seqs = brca1_aln,
    group = "Z64",
    cube = c("ACGT", "TGCA"),
    cube_alt = c("CATG", "GTAC"),
    nms = nams
)
usethis::use_data(brca1_autm, brca1_aln, overwrite = TRUE, compress = "xz")

autby_coef <- automorphism_bycoef(brca1_autm)
autby_coef

usethis::use_data(autby_coef, overwrite = TRUE, compress = "xz")


nams <- c(paste0("human_1.", 0:21),"human_2","gorilla_1","gorilla_2","gorilla_3",
          "chimpanzee_1","chimpanzee_2","chimpanzee_3","chimpanzee_4",
          "bonobos_1","bonobos_2","bonobos_3","bonobos_4","silvery_gibbon_1",
          "silvery_gibbon_1","silvery_gibbon_3","golden_monkey_1",
          "golden_monkey_2","gelada_baboon","bolivian_monkey")

URL <- paste0("https://github.com/genomaths/seqalignments/raw/master/BRCA1/",
              "brca1_primates_dna_repair_41_sequences.fasta")

brca1_aln2 <- readDNAMultipleAlignment(filepath = URL)

brca1_autm2 <- automorphisms(
    seqs = brca1_aln2,
    group = "Z64",
    cube = c("ACGT", "TGCA"),
    cube_alt = c("CATG", "GTAC"),
    nms = nams
)
usethis::use_data(brca1_autm2, brca1_aln2, overwrite = TRUE, compress = "xz")



URL <- paste0(
    "https://github.com/genomaths/seqalignments/raw/master/",
    "COVID-19/AY390556.1_265-13398_13398-21485_RNA-POL_SARS_COVI_GZ02.fas"
)

autm <- automorphisms(
    filepath = URL,
    group = "Z64",
    cube = c("ACGT", "TGCA"),
    cube_alt = c("CATG", "GTAC")
)

usethis::use_data(autm, overwrite = TRUE)


URL <- paste0(
    "https://github.com/genomaths/seqalignments/raw/master/",
    "COVID-19/AY390556.1_and_KY417151.1_aligned_protein-coding.fas"
)

URL <- paste0(
    "https://github.com/genomaths/seqalignments/raw/master/", 
    "COVID-19/AY390556.1_and_KY417151.1_aligned_protein-coding.fas")

covid_aln <- readDNAMultipleAlignment(filepath = URL)
covid_aln

usethis::use_data(covid_aln, overwrite = TRUE)


covid_autm <- automorphisms(
    seq = covid_aln,
    group = "Z64",
    cube = c("ACGT", "TGCA"),
    cube_alt = c("CATG", "GTAC")
)
covid_autm

usethis::use_data(covid_autm, overwrite = TRUE)

autm_z125 <- automorphisms(
    seq = covid_aln, 
    group = "Z125", 
    cube = c("ACGT", "TGCA"),
    cube_alt = c("CATG", "GTAC"),
    verbose = FALSE)
autm_z125
usethis::use_data(autm_z125, overwrite = TRUE)

autm_3d <- automorphisms(
    seq = covid_aln, 
    group = "Z5^3", 
    cube = c("ACGT", "TGCA"),
    cube_alt = c("CATG", "GTAC"),
    verbose = FALSE)
autm_3d
usethis::use_data(autm_3d, overwrite = TRUE)


## Codon distance matrices for the standard genetic code on Z4
cube = c("ACGT", "AGCT", "TCGA", "TGCA", "CATG", 
         "GTAC", "CTAG", "GATC", "ACTG", "ATCG", 
         "GTCA", "GCTA", "CAGT", "TAGC", "TGAC", 
         "CGAT", "AGTC", "ATGC", "CGTA", "CTGA", 
         "GACT", "GCAT", "TACG", "TCAG")

names(cube) <- cube

cdm_z64 <- lapply(cube, function(x) 
                codon_dist_matrix(
                    cube = x,
                    group = "Z4", 
                    output = "vector",
                    num.cores = 20L))

usethis::use_data(cdm_z64, overwrite = TRUE, compress = "xz")


