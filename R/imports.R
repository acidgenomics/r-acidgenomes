## Classes =====================================================================

#' @importClassesFrom GenomicRanges GenomicRanges GenomicRangesList
#' @importClassesFrom IRanges DataFrameList
#' @importClassesFrom S4Vectors DFrame DataFrame Vector
NULL



## S4 generics and methods =====================================================

#' @importFrom AcidGenerics as.DataFrame camelCase encode kebabCase leftJoin
#' makeNames mutateAt rbindToDataFrame removeNA sanitizeNA snakeCase
#' upperCamelCase
#' @importFrom BiocGenerics %in% append as.data.frame do.call grep grepl
#' is.unsorted lapply match order rbind setdiff sort unlist
#' @importFrom GenomeInfoDb genome genome<- seqinfo seqinfo<- seqlevels
#' @importFrom IRanges gsub ranges trim
#' @importFrom S4Vectors complete.cases decode expand head mcols mcols<-
#' metadata metadata<- na.omit split tail
#' @importFrom pipette export import
#'
#' @importMethodsFrom AcidPlyr leftJoin mutateAt rbindToDataFrame
#' @importMethodsFrom pipette as.DataFrame export import removeNA sanitizeNA
#' @importMethodsFrom syntactic camelCase kebabCase makeNames snakeCase
#' upperCamelCase
NULL



## Standard functions ==========================================================

#' @importFrom AcidBase download fileExt forceDetach initDir parentDir pasteURL
#' printString realpath requireNamespaces showSlotInfo standardizeCall
#' @importFrom AcidCLI abort alert alertInfo alertSuccess alertWarning dl h1
#' toInlineString
#' @importFrom GenomeInfoDb Seqinfo getChromInfoFromEnsembl
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges CharacterList DataFrameList IntegerList IRanges
#' SplitDataFrameList
#' @importFrom S4Vectors DataFrame Rle
#' @importFrom goalie allAreAtomic allAreFiles allAreMatchingFixed
#' allAreMatchingRegex allAreNotMatchingRegex allAreURLs areDisjointSets
#' areIntersectingSets areSetEqual assert bapply hasColnames hasDuplicates
#' hasInternet hasLength hasNoDuplicates hasNames hasRownames hasRows
#' hasValidNames isADir isAFile isAURL isAny isCharacter isFlag isInstalled
#' isInt isMatchingFixed isMatchingRegex isNotMatchingFixed isOrganism
#' isScalar isString isSubset isSystemCommand isWindows validate
#' validateClasses
#' @importFrom methods as is isClass new setClass setGeneric setMethod
#' setValidity signature validObject
#' @importFrom pipette cacheURL getJSON getURLDirList md5 sha256
#' @importFrom stringi stri_match_first_regex
#' @importFrom utils capture.output packageName packageVersion sessionInfo
NULL
