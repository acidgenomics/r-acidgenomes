## FIXME NEED TO REWORK IDENTIFIER MATCHING FOR SYNONYMS...
## CURRENTLY FAILING IN DEPMAPANALYSIS DRAFT UPDATE.



#' AcidGenomes
#'
#' Toolkit for downloading and processing genome annotations.
#'
#' @keywords internal
#'
#' @importClassesFrom GenomicRanges CompressedGRangesList GRanges GRangesList
#' @importClassesFrom IRanges DataFrameList
#' @importClassesFrom S4Vectors DataFrame Vector
#'
#' @importFrom AcidBase compress download fileExt forceDetach initDir pasteURL
#'   printString realpath requireNamespaces showSlotInfo
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
#'   seqinfo<- seqlengths seqlevels seqnames
#' @importFrom GenomicRanges GRanges ranges trim
#' @importFrom IRanges DataFrameList IRanges
#' @importFrom S4Vectors DataFrame Rle
#' @importFrom ensembldb ensemblVersion genes listColumns transcripts
#' @importFrom goalie allAreAtomic allAreFiles allAreMatchingRegex
#'   allAreNotMatchingRegex allAreURLs areDisjointSets areIntersectingSets
#'   areSetEqual assert bapply hasColnames hasDuplicates hasInternet hasLength
#'   hasNoDuplicates hasNames hasRownames hasRows hasValidNames isADir isAFile
#'   isAURL isAny isCharacter isFlag isInstalled isInt isMatchingFixed
#'   isMatchingRegex isNotMatchingFixed isOrganism isScalar isString isSubset
#'   isSystemCommand isWindows validate validateClasses
#' @importFrom httr GET content content_type stop_for_status
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom methods as is isClass new setClass setGeneric setMethod
#'   setValidity signature validObject
#' @importFrom pipette as_tibble cacheURL export getJSON getURL getURLDirList
#'   import md5 removeNA sanitizeNA sha256
#' @importFrom purrr map_df
#' @importFrom stringr boundary str_extract str_match str_split_fixed
#' @importFrom syntactic camelCase kebabCase makeNames snakeCase upperCamelCase
#' @importFrom utils capture.output packageName packageVersion sessionInfo
#'
#' @importMethodsFrom GenomicRanges is.unsorted sort
"_PACKAGE"
