for (organism in c(
    ## > Arabidopsis thaliana",
    "Bos taurus",
    "Caenorhabditis elegans",
    "Canis lupus familiaris",
    "Danio rerio",
    "Drosophila melanogaster",
    "Gallus gallus",
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Saccharomyces cerevisiae",
    "Sus scrofa"
)) {
    test_that(organism, {
        object <- goTermsPerGeneName(organism)
        expect_s4_class(object, "DFrame")
    })
}
