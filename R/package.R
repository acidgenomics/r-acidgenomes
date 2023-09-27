#' AcidGenomes
#'
#' Toolkit for downloading and processing genome annotations.
#'
#' @keywords internal
"_PACKAGE"



## Classes =====================================================================

#' @importClassesFrom GenomicRanges CompressedGRangesList GRanges GRangesList
#' @importClassesFrom IRanges DFrameList
#' @importClassesFrom S4Vectors DFrame Vector
NULL



## S4 generics and methods =====================================================

#' @importFrom AcidGenerics EnsemblToNcbi GeneToSymbol NcbiToEnsembl TxToGene
#' as.DataFrame camelCase encode export import kebabCase leftJoin makeNames
#' matchNested mutateAt rbindToDataFrame removeNA sanitizeNA snakeCase
#' stripGeneVersions stripTranscriptVersions upperCamelCase
#' @importFrom BiocGenerics %in% append as.data.frame do.call grep grepl
#' is.unsorted lapply match order organism organism<- rbind setdiff sort unlist
#' updateObject
#' @importFrom GenomeInfoDb genome genome<- seqinfo seqinfo<- seqlevels
#' seqnames<- seqnames
#' @importFrom IRanges gsub ranges trim
#' @importFrom S4Vectors complete.cases decode expand head mcols mcols<-
#' merge metadata metadata<- na.omit split summary tail
NULL

#' @importMethodsFrom AcidBase matchNested
#' @importMethodsFrom AcidPlyr leftJoin mutateAt rbindToDataFrame
#' @importMethodsFrom pipette as.DataFrame export import removeNA sanitizeNA
#' @importMethodsFrom syntactic camelCase kebabCase makeNames snakeCase
#' upperCamelCase
NULL



## Standard functions ==========================================================

#' @importFrom AcidBase download dupes fileExt initDir parentDir pasteURL
#' quietly realpath showHeader showSlotInfo standardizeCall strExtract strMatch
#' strSplit
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
#' hasValidNames isADir isAFile isAURL isAnExistingURL isAny isCharacter
#' isDuplicate isFlag isInstalled isInt isIntegerish isMatchingFixed
#' isMatchingRegex isOrganism isScalar isString isSubset isSystemCommand
#' isWindows requireNamespaces validate validateClasses
#' @importFrom methods as is isClass new setClass setGeneric setMethod
#' setValidity show signature validObject
#' @importFrom pipette cacheURL fillLines getJSON getURLDirList md5 sha256
#' @importFrom utils capture.output packageName packageVersion sessionInfo
NULL
