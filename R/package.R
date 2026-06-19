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
#' @importFrom AcidGenerics as.DataFrame camelCase export import kebabCase
#' @importFrom AcidGenerics leftJoin makeNames matchNested mutateAt
#' @importFrom AcidGenerics rbindToDataFrame removeNa sanitizeNa snakeCase
#' @importFrom AcidGenerics stripExonVersions stripGeneVersions
#' @importFrom AcidGenerics stripTranscriptVersions upperCamelCase
#' @importFrom BiocGenerics %in% append as.data.frame do.call grep grepl
#' @importFrom BiocGenerics is.unsorted lapply match order organism
#' @importFrom BiocGenerics organism<- rbind setdiff sort unlist updateObject
#' @importFrom GenomeInfoDb genome genome<- seqinfo seqinfo<- seqlevels
#' @importFrom GenomeInfoDb seqnames<- seqnames
#' @importFrom IRanges gsub ranges startsWith sub toupper trim
#' @importFrom S4Vectors complete.cases decode expand head mcols mcols<-
#' @importFrom S4Vectors merge metadata metadata<- na.omit split summary tail
NULL

#' @importMethodsFrom AcidBase matchNested
#' @importMethodsFrom AcidPlyr leftJoin mutateAt rbindToDataFrame
#' @importMethodsFrom pipette as.DataFrame export import removeNa sanitizeNa
#' @importMethodsFrom syntactic camelCase kebabCase makeNames snakeCase
#' @importMethodsFrom syntactic upperCamelCase
NULL


## Standard functions ==========================================================

#' @importFrom AcidBase download dupes fileExt initDir matchAll parentDir
#' @importFrom AcidBase pasteUrl quietly realpath showHeader showSlotInfo
#' @importFrom AcidBase standardizeCall strExtract strMatch strSplit
#' @importFrom AcidCLI abort alert alertInfo alertSuccess alertWarning dl h1
#' @importFrom AcidCLI toInlineString
#' @importFrom GenomeInfoDb Seqinfo getChromInfoFromEnsembl
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges CharacterList DataFrameList IntegerList IRanges
#' @importFrom IRanges SplitDataFrameList
#' @importFrom S4Vectors DataFrame
#' @importFrom goalie allAreAtomic allAreFiles allAreMatchingFixed
#' @importFrom goalie allAreMatchingRegex allAreNotMatchingRegex allAreUrls
#' @importFrom goalie areDisjointSets areIntersectingSets areSetEqual assert
#' @importFrom goalie bapply hasColnames hasDuplicates hasInternet hasLength
#' @importFrom goalie hasNoDuplicates hasNames hasRownames hasRows hasValidNames
#' @importFrom goalie isADir isAFile isAUrl isAnExistingUrl isAny isCharacter
#' @importFrom goalie isDuplicate isFlag isInstalled isInt isIntegerish
#' @importFrom goalie isMatchingFixed isMatchingRegex isOrganism isScalar
#' @importFrom goalie isString isSubset isSystemCommand isWindows
#' @importFrom goalie requireNamespaces validate validateClasses
#' @importFrom methods as is isClass new setClass setGeneric setMethod
#' @importFrom methods setValidity show signature validObject
#' @importFrom parallel mclapply
#' @importFrom pipette cacheUrl getJson getUrlDirList md5 sha256
#' @importFrom utils capture.output packageName packageVersion sessionInfo
#' @importFrom utils write.table
NULL
