#' AcidGenomes
#'
#' Toolkit for downloading and processing genome annotations.
#'
#' @keywords internal
#'
#' @importClassesFrom AcidGenerics DataFrame DataFrameList GenomicRanges
#'   GenomicRangesList Vector
#'
#' @importMethodsFrom pipette as.DataFrame coerce
#'
#' @importFrom AcidBase capture.output download fileExt forceDetach initDir
#'   packageName packageVersion parentDir pasteURL printString realpath
#'   requireNamespaces sessionInfo showSlotInfo standardizeCall
#' @importFrom AcidCLI abort alert alertInfo alertSuccess alertWarning dl h1
#'   toInlineString
#' @importFrom AcidGenerics CharacterList DataFrame DataFrameList GRanges
#'   IntegerList IRanges Rle Seqinfo SplitDataFrameList %in% append as.DataFrame
#'   as.data.frame complete.cases decode do.call encode expand genome genome<-
#'   getChromInfoFromEnsembl grep grepl gsub head is.unsorted lapply match mcols
#'   mcols<- metadata metadata<- na.omit order ranges rbind setdiff seqinfo
#'   seqinfo<- seqlevels sort split tail trim unlist
#' @importFrom AcidPlyr leftJoin mutateAt rbindToDataFrame
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
#' @importFrom pipette as_tibble cacheURL export getJSON getURL getURLDirList
#'   import md5 removeNA sanitizeNA sha256
#' @importFrom stringr boundary str_extract str_match str_split_fixed
#' @importFrom syntactic camelCase kebabCase makeNames snakeCase upperCamelCase
"_PACKAGE"



#' @name params
#' @inherit AcidRoxygen::params return title
#' @keywords internal
#'
#' @param cache `logical(1)`.
#'   Cache URLs locally, using BiocFileCache internally.
#' @param ignoreVersion `logical(1)`.
#'   Ignore identifier (e.g. transcript, gene) versions.
#'   When applicable, the identifier containing version numbers will be stored
#'   in `txIdVersion` and `geneIdVersion`, and the variants without versions
#'   will be stored in `txId`, `txIdNoVersion`, `geneId`, and `geneIdNoVersion`.
#' @param synonyms `logical(1)`.
#'   Include gene synonyms.
#'   Queries the Ensembl web server, and is CPU intensive.
NULL
