#' AcidGenomes
#'
#' Toolkit for downloading and processing genome annotations.
#'
#' @keywords internal
#'
#' @importClassesFrom AcidGenerics DataFrame DataFrameList Vector
#' @importClassesFrom pipette CompressedGRangesList GRangesList GRanges
#'
#' @importFrom AcidBase capture.output download fileExt forceDetach initDir
#'   packageName packageVersion parentDir pasteURL printString realpath
#'   requireNamespaces sessionInfo showSlotInfo standardizeCall
#' @importFrom AcidCLI abort alert alertInfo alertSuccess alertWarning dl h1
#'   toInlineString
#' @importFrom AcidGenerics CharacterList DataFrame DataFrameList IntegerList
#'   IRanges Rle SplitDataFrameList %in% append as.data.frame complete.cases
#'   decode do.call encode expand grep grepl gsub head is.unsorted lapply match
#'   mcols mcols<- metadata metadata<- na.omit order ranges rbind setdiff sort
#'   split tail trim unlist
#' @importFrom AcidPlyr leftJoin mutateAt rbindToDataFrame
#' @importFrom GenomeInfoDb Seqinfo getChromInfoFromEnsembl genome genome<-
#'   seqinfo seqinfo<- seqlevels
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
#' @importFrom pipette GRanges as.DataFrame as_tibble cacheURL export getJSON
#'   getURL getURLDirList import md5 removeNA sanitizeNA sha256
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
