## Classes =====================================================================

#' @importClassesFrom GenomicRanges GenomicRanges GenomicRangesList
#' @importClassesFrom IRanges DataFrameList
#' @importClassesFrom S4Vectors DFrame DataFrame Vector
NULL



## S4 generics and methods =====================================================

#' @importFrom AcidGenerics Ensembl2Ncbi Gene2Symbol Ncbi2Ensembl Tx2Gene
#' as.DataFrame camelCase encode kebabCase leftJoin makeNames matchNested
#' mutateAt rbindToDataFrame removeNA sanitizeNA snakeCase stripGeneVersions
#' stripTranscriptVersions upperCamelCase
#' @importFrom BiocGenerics %in% append as.data.frame do.call grep grepl
#' is.unsorted lapply match order organism organism<- rbind setdiff sort unlist
#' updateObject
#' @importFrom GenomeInfoDb genome genome<- seqinfo seqinfo<- seqlevels
#' seqnames<- seqnames
#' @importFrom IRanges gsub ranges trim
#' @importFrom S4Vectors complete.cases decode expand head mcols mcols<-
#' metadata metadata<- na.omit split summary tail
#' @importFrom pipette export import
NULL

#' @importMethodsFrom AcidBase matchNested
#' @importMethodsFrom AcidPlyr leftJoin mutateAt rbindToDataFrame
#' @importMethodsFrom pipette as.DataFrame export import removeNA sanitizeNA
#' @importMethodsFrom syntactic camelCase kebabCase makeNames snakeCase
#' upperCamelCase
NULL



## Standard functions ==========================================================

#' @importFrom AcidBase download fileExt forceDetach initDir parentDir pasteURL
#' printString realpath showSlotInfo standardizeCall
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
#' isInt isIntegerish isMatchingFixed isMatchingRegex isOrganism isScalar
#' isString isSubset isSystemCommand isWindows requireNamespaces validate
#' validateClasses
#' @importFrom methods as is isClass new setClass setGeneric setMethod
#' setValidity signature validObject
#' @importFrom pipette cacheURL getJSON getURLDirList md5 sha256
#' @importFrom stringi stri_match_first_regex
#' @importFrom utils capture.output packageName packageVersion sessionInfo
NULL
