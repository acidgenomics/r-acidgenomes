#' All global variables
#' @noRd
NULL



## Updated 2021-01-22.
.gffFormats <- c("GFF3", "GTF")



## Updated 2021-01-22.
.rtracklayerFormats <- c("GFF3", "GTF")



## Updated 2021-01-22.
.rtracklayerProviders <-
    c("Ensembl", "FlyBase", "GENCODE", "RefSeq", "UCSC", "WormBase")



#' Package version
#'
#' @note Updated 2020-10-06.
#' @noRd
.version <- packageVersion(packageName())



#' AcidGenomes test data URL
#'
#' @export
#' @keywords internal
#' @note Updated 2020-10-06.
#'
#' @examples
#' AcidGenomesTestsURL
AcidGenomesTestsURL <-  # nolint
    paste0(
        "https://tests.acidgenomics.com/AcidGenomes/",
        "v", .version$major, ".", .version$minor  # nolint
    )
