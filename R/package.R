#' AcidGenomes
#'
#' Toolkit for downloading and processing genome annotations.
#'
#' @keywords internal
#'
#' @importClassesFrom AcidGenerics DataFrame DataFrameList Vector
#' @importClassesFrom pipette CompressedGRangesList GRangesList GRanges
#'
#' @importFrom AcidBase capture.output compress download fileExt forceDetach
#'   initDir packageName packageVersion pasteURL printString realpath
#'   requireNamespaces sessionInfo showSlotInfo standardizeCall
#' @importFrom AcidCLI alert alertInfo alertSuccess alertWarning dl h1
#'   toInlineString
#' @importFrom AcidGenerics CharacterList DataFrame DataFrameList IntegerList
#'   IRanges Rle %in% append as.data.frame complete.cases decode do.call encode
#'   expand grep grepl gsub head is.unsorted lapply match mcols mcols<- metadata
#'   metadata<- na.omit order ranges setdiff sort split tail trim unlist
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
#' @importFrom pipette GRanges as_tibble cacheURL export getJSON getURL
#'   getURLDirList import md5 removeNA sanitizeNA sha256
#' @importFrom stringr boundary str_extract str_match str_split_fixed
#' @importFrom syntactic camelCase kebabCase makeNames snakeCase upperCamelCase
"_PACKAGE"
