#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges DataFrameList
#' @importClassesFrom S4Vectors DataFrame Vector
#'
#' @importFrom AcidBase forceDetach matchArgsToDoCall pasteURL printString
#'   requireNamespaces
#' @importFrom AcidPlyr leftJoin
#' @importFrom AnnotationDbi select
#' @importFrom AnnotationHub AnnotationHub query snapshotDate
#' @importFrom BiocParallel bplapply
#' @importFrom GenomicRanges GRanges ranges seqnames
#' @importFrom IRanges DataFrameList IRanges
#' @importFrom S4Vectors DataFrame Rle %in% as.data.frame complete.cases decode
#'   do.call expand head lapply match mcols mcols<- metadata metadata<- order
#'   setdiff split summary tail
#' @importFrom SummarizedExperiment rowData rowRanges
#' @importFrom cli cli_alert cli_alert_info cli_alert_warning
#'   cli_div cli_dl cli_end
#' @importFrom ensembldb ensemblVersion genes transcripts
#' @importFrom goalie areDisjointSets areSetEqual assert bapply hasDuplicates
#'   hasInternet hasLength hasNoDuplicates hasNames hasRownames hasRows isAny
#'   isCharacter isFlag isInt isMatchingRegex isString isSubset validate
#' @importFrom methods as is new setClass setGeneric setMethod setValidity
#'   signature
#' @importFrom pipette as_tibble as.SummarizedExperiment cacheURL getURL
#'   getURLDirList import removeNA sanitizeNA
#' @importFrom stringr boundary str_extract str_match str_split_fixed
#' @importFrom syntactic camelCase makeNames
#' @importFrom utils capture.output packageName packageVersion
NULL
