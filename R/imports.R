#' @importClassesFrom GenomicRanges CompressedGRangesList GRanges
#' @importClassesFrom IRanges DataFrameList
#' @importClassesFrom S4Vectors DataFrame Vector
#'
#' @importFrom AcidBase compress download fileExt initDir matchArgsToDoCall
#'   pasteURL printString realpath requireNamespaces showSlotInfo
#' @importFrom AcidCLI alert alertInfo alertSuccess alertWarning dl h1
#'   toInlineString
#' @importFrom AcidGenerics %in% as.data.frame complete.cases decode do.call
#'   expand head is.unsorted lapply match mcols mcols<- metadata metadata<-
#'   order setdiff sort split tail
#' @importFrom AcidPlyr leftJoin
#' @importFrom AnnotationDbi columns select
#' @importFrom AnnotationHub AnnotationHub query snapshotDate
#' @importFrom BiocParallel bplapply
#' @importFrom GenomeInfoDb Seqinfo getChromInfoFromEnsembl genome isCircular
#'   seqinfo<- seqlengths seqnames
#' @importFrom GenomicRanges GRanges ranges trim
#' @importFrom IRanges DataFrameList IRanges
#' @importFrom S4Vectors DataFrame Rle
#' @importFrom SummarizedExperiment rowData rowRanges
#' @importFrom digest digest
#' @importFrom ensembldb ensemblVersion genes listColumns transcripts
#' @importFrom goalie allAreFiles allAreMatchingRegex allAreNotMatchingRegex
#'   allAreURLs areDisjointSets areIntersectingSets areSetEqual assert bapply
#'   hasDuplicates hasInternet hasLength hasNoDuplicates hasNames hasRownames
#'   hasRows isADir isAFile isAURL isAny isCharacter isFlag isInt
#'   isMatchingFixed isMatchingRegex isNotMatchingFixed isOrganism isScalar
#'   isString isSubset isSystemCommand isWindows validate
#' @importFrom httr GET content content_type stop_for_status
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom methods as is isClass new setClass setGeneric setMethod
#'   setValidity signature
#' @importFrom pipette as_tibble as.SummarizedExperiment cacheURL export
#'   getJSON getURL getURLDirList import removeNA sanitizeNA
#' @importFrom purrr map_df
#' @importFrom stringr boundary str_extract str_match str_split_fixed
#' @importFrom syntactic camelCase kebabCase makeNames snakeCase upperCamelCase
#' @importFrom utils capture.output packageName packageVersion sessionInfo
#'
#' @importMethodsFrom GenomicRanges is.unsorted sort
NULL
