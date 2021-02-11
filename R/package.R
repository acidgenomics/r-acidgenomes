## FIXME NEED TO REWORK IDENTIFIER MATCHING FOR SYNONYMS...
## CURRENTLY FAILING IN DEPMAPANALYSIS DRAFT UPDATE.

## FIXME CAN WE TAKE BIOCPARALLEL OUT COMPLETELY HERE?



#' AcidGenomes
#'
#' Toolkit for downloading and processing genome annotations.
#'
#' @keywords internal
#'
#' @importClassesFrom AcidGenerics DataFrame Vector
#' @importClassesFrom pipette DataFrameList
#'
#' @importFrom AcidBase capture.output compress download fileExt forceDetach
#'   initDir packageName packageVersion pasteURL printString realpath
#'   requireNamespaces sessionInfo showSlotInfo
#' @importFrom AcidCLI alert alertInfo alertSuccess alertWarning dl h1
#'   toInlineString
#' @importFrom AcidGenerics DataFrame Rle %in% as.data.frame complete.cases
#'   decode do.call expand head is.unsorted lapply match mcols mcols<- metadata
#'   metadata<- order setdiff sort split tail
#' @importFrom AcidPlyr leftJoin map_df
#' @importFrom AnnotationDbi columns select
#' @importFrom AnnotationHub AnnotationHub query snapshotDate
#' @importFrom GenomeInfoDb Seqinfo getChromInfoFromEnsembl genome isCircular
#'   seqinfo<- seqlengths seqlevels seqnames
#' @importFrom ensembldb ensemblVersion genes listColumns transcripts
#' @importFrom goalie allAreAtomic allAreFiles allAreMatchingFixed
#'   allAreMatchingRegex allAreNotMatchingRegex allAreURLs areDisjointSets
#'   areIntersectingSets areSetEqual assert bapply hasColnames hasDuplicates
#'   hasInternet hasLength hasNoDuplicates hasNames hasRownames hasRows
#'   hasValidNames isADir isAFile isAURL isAny isCharacter isFlag isInstalled
#'   isInt isMatchingFixed isMatchingRegex isNotMatchingFixed isOrganism
#'   isScalar isString isSubset isSystemCommand isWindows validate
#'   validateClasses
#' @importFrom httr GET content content_type stop_for_status

#' @importFrom methods as is isClass new setClass setGeneric setMethod
#'   setValidity signature validObject
#' @importFrom pipette DataFrameList GRanges GRangesList IRanges as_tibble
#'   cacheURL export getJSON getURL getURLDirList import md5 ranges removeNA
#'   sanitizeNA sha256 trim
#' @importFrom stringr boundary str_extract str_match str_split_fixed
#' @importFrom syntactic camelCase kebabCase makeNames snakeCase upperCamelCase
"_PACKAGE"
