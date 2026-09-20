test_that("Ensembl data URLs use HTTPS", {
    expect_identical(
        object = .ensemblFtpUrl("release-116", "gtf"),
        expected = "https://ftp.ensembl.org/pub/release-116/gtf"
    )
})

test_that("Ensembl release versions are parsed strictly", {
    expect_identical(
        object = .parseEnsemblVersion(" 116 "),
        expected = 116L
    )
    expect_error(
        object = .parseEnsemblVersion("Ensembl Release 116"),
        regexp = "Failed to extract release version"
    )
})
