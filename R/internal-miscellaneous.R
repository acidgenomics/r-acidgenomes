#' Download multiple genome files in a single call
#'
#' @note Updated 2022-05-24.
#' @noRd
#'
#' @return `character`
#' Local file paths.
.downloadUrls <-
    function(urls, outputDir, cache) {
        assert(
            allAreUrls(urls),
            isString(outputDir),
            isFlag(cache)
        )
        outputDir <- initDir(outputDir)
        destFiles <- vapply(
            X = urls,
            FUN = function(url) {
                file.path(outputDir, basename(url))
            },
            FUN.VALUE = character(1L),
            USE.NAMES = TRUE
        )
        assert(identical(names(urls), names(destFiles)))
        if (isTRUE(cache)) {
            cacheFiles <- .cacheIt(file = urls)
            file.copy(
                from = cacheFiles,
                to = destFiles,
                overwrite = TRUE,
                recursive = FALSE
            )
        } else {
            Map(url = urls, destfile = destFiles, f = download)
        }
        out <- destFiles
        names(out) <- names(urls)
        assert(allAreFiles(out))
        invisible(out)
    }


#' Get an internal function from the package
#'
#' @note Updated 2021-03-03.
#' @noRd
.getFun <- function(x) {
    fun <- get(
        x = x,
        envir = asNamespace(.pkgName),
        inherits = FALSE
    )
    assert(is.function(fun))
    fun
}


#' Map a genome build to UCSC
#'
#' @details
#' Currently used internally for handoff to `Seqinfo` (e.g. GENCODE genome).
#'
#' @note Updated 2022-05-04.
#' @noRd
#'
#' @examples
#' .mapGenomeBuildToUCSC("GRCh38")
.mapGenomeBuildToUCSC <- function(x) {
    assert(isString(x))
    switch(
        EXPR = x,
        BDGP5 = "dm3", # Drosophila melanogaster
        BDGP6 = "dm6", # Drosophila melanogaster
        CanFam3.1 = "canFam3", # Canis familiaris
        GRCg6a = "galGal6", # Gallus gallus
        GRCh37 = "hg19", # Homo sapiens
        GRCh38 = "hg38", # Homo sapiens
        GRCm37 = "mm9", # Mus musculus
        GRCm38 = "mm10", # Mus musculus
        GRCm39 = "mm39", # Mus musculus
        GRCz10 = "danRer10", # Danio rerio
        GRCz11 = "danRer11", # Danio rerio
        Galgal4 = "galGal4", # Gallus gallus
        "Gallus_gallus-5.0" = "galGal5", # Gallus gallus
        JGI_4.1 = "xenTro2", # Xenopus tropicalis
        JGI_4.2 = "xenTro3", # Xenopus tropicalis
        "R64-1-1" = "sacCer3", # Saccharomyces cerevisiae
        Rnor_5.0 = "rn5", # Rattus norvegicus
        Rnor_6.0 = "rn6", # Rattus norvegicus
        Sscrofa11.1 = "susScr11", # Sus scrofa
        TAIR10 = "araTha1", # Arabidopsis thaliana
        WBcel235 = "ce11", # Caenorhabditis elegans
        WS220 = "ce10", # Caenorhabditis elegans
        abort(sprintf("Unsupported genome build: {.val %s}.", x))
    )
}


#' Nest an S4 Data Frame by a grouping factor
#'
#' Consider migrating this to AcidPlyr in a future update.
#'
#' @note Updated 2023-04-26.
#' @noRd
#'
#' @param object Object.
#' @param by `character(1)`.
#' Identifier column.
#' @param exclude `character` or `NULL`.
#' Column names to exclude.
#'
#' @return `DFrame`.
#'
#' @seealso
#' - `tidyr::nest`.
.nest2 <- function(object, by, exclude = NULL) {
    assert(
        isString(by),
        isCharacter(exclude, nullOk = TRUE)
    )
    if (!is.null(exclude)) {
        object <- object[, setdiff(colnames(object), exclude)]
    }
    object <- unique(object[complete.cases(object), ])
    spl <- split(x = object, f = as.factor(object[[by]]))
    args <- list()
    args[[by]] <- names(spl)
    args <- append(
        x = args,
        values = Map(
            col = setdiff(colnames(object), by),
            MoreArgs = list(spl = spl),
            f = function(col, spl) {
                unname(spl[, col])
            }
        )
    )
    do.call(what = DataFrame, args = args)
}


#' Generate chromosome sizes file from a genome FASTA
#'
#' Produces a 2-column tab-delimited file (seqname, length) suitable for use
#' with kallisto index `--genomebam` / `--chromosomes` and similar tools.
#' Parses FASTA headers from the compressed file using only base R connections.
#'
#' @note Updated 2026-05-31.
#' @noRd
#'
#' @param fastaFile `character(1)`.
#' Path to a genome FASTA file (`.fa.gz` or `.fa`).
#'
#' @param outputDir `character(1)`.
#' Output directory where `chrom.sizes` will be written.
#'
#' @return Invisible `character(1)`. Path to the written `chrom.sizes` file.
.generateChromSizes <- function(fastaFile, outputDir) {
    assert(
        isAFile(fastaFile),
        isString(outputDir)
    )
    alert("Generating chromosome sizes file for kallisto.")
    con <- if (endsWith(fastaFile, ".gz")) {
        gzcon(file(fastaFile, open = "rb"))
    } else {
        file(fastaFile, open = "r")
    }
    on.exit(close(con), add = TRUE)
    seqNames <- character(0L)
    seqLengths <- integer(0L)
    currentSeq <- NULL
    currentLen <- 0L
    repeat {
        line <- readLines(con = con, n = 1L, warn = FALSE)
        if (length(line) == 0L) {
            break
        }
        if (startsWith(line, ">")) {
            if (!is.null(currentSeq)) {
                seqNames <- c(seqNames, currentSeq)
                seqLengths <- c(seqLengths, currentLen)
            }
            ## Take only the first token from the FASTA header.
            currentSeq <- sub("^>([^ \t]+).*$", "\\1", line)
            currentLen <- 0L
        } else {
            currentLen <- currentLen + nchar(line)
        }
    }
    if (!is.null(currentSeq)) {
        seqNames <- c(seqNames, currentSeq)
        seqLengths <- c(seqLengths, currentLen)
    }
    assert(
        hasLength(seqNames),
        hasLength(seqLengths),
        identical(length(seqNames), length(seqLengths))
    )
    outputFile <- file.path(outputDir, "chrom.sizes")
    df <- data.frame(
        seqname = seqNames,
        length = seqLengths,
        stringsAsFactors = FALSE
    )
    write.table(
        x = df,
        file = outputFile,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
    alertSuccess(sprintf(
        "Chromosome sizes file written to {.path %s}.",
        outputFile
    ))
    invisible(outputFile)
}


#' Generate salmon decoy files (gentrome + decoys.txt)
#'
#' Salmon selective alignment recommends indexing a "gentrome" — the
#' transcriptome FASTA concatenated with the genome FASTA — along with a
#' `decoys.txt` file listing the genome sequence names. This avoids mapping
#' reads to non-transcriptomic regions.
#'
#' @note Updated 2026-05-31.
#' @noRd
#'
#' @param genomeFasta `character(1)`.
#' Path to the genome FASTA (`.fa.gz` or `.fa`).
#'
#' @param transcriptomeFasta `character(1)`.
#' Path to the transcriptome FASTA (`.fa.gz` or `.fa`).
#'
#' @param outputDir `character(1)`.
#' Output directory where `gentrome.fa.gz` and `decoys.txt` will be written.
#'
#' @return Invisible named `character` with paths to `gentrome` and `decoys`.
.generateSalmonDecoys <- function(genomeFasta, transcriptomeFasta, outputDir) {
    assert(
        isAFile(genomeFasta),
        isAFile(transcriptomeFasta),
        isString(outputDir)
    )
    alert("Generating salmon decoy files (gentrome + decoys.txt).")
    outputDir <- initDir(outputDir)
    ## Extract genome sequence names for decoys.txt.
    con <- if (endsWith(genomeFasta, ".gz")) {
        gzcon(file(genomeFasta, open = "rb"))
    } else {
        file(genomeFasta, open = "r")
    }
    on.exit(close(con), add = TRUE)
    decoyNames <- character(0L)
    repeat {
        line <- readLines(con = con, n = 1L, warn = FALSE)
        if (length(line) == 0L) {
            break
        }
        if (startsWith(line, ">")) {
            decoyNames <- c(decoyNames, sub("^>([^ \t]+).*$", "\\1", line))
        }
    }
    close(con)
    on.exit(NULL)
    assert(hasLength(decoyNames))
    decoyFile <- file.path(outputDir, "decoys.txt")
    writeLines(text = decoyNames, con = decoyFile)
    ## Concatenate transcriptome + genome into gentrome.fa.gz using system tools
    ## when available, otherwise fall back to R connections for portability.
    gentromeFile <- file.path(outputDir, "gentrome.fa.gz")
    if (!isWindows() && nzchar(Sys.which("cat")) && nzchar(Sys.which("gzip"))) {
        ## Fast path: shell pipeline (handles both .gz and plain FASTA inputs).
        decompressCmd <- function(f) {
            if (endsWith(f, ".gz")) {
                paste("gzip -dc", shQuote(f))
            } else {
                paste("cat", shQuote(f))
            }
        }
        cmd <- paste(
            decompressCmd(transcriptomeFasta),
            decompressCmd(genomeFasta),
            "| gzip -c >",
            shQuote(gentromeFile)
        )
        ret <- system(cmd)
        if (!identical(ret, 0L)) {
            abort("Failed to generate gentrome.fa.gz via shell pipeline.")
        }
    } else {
        ## Portable R fallback: read both FASTAs and write gzipped output.
        outCon <- gzfile(gentromeFile, open = "wb")
        on.exit(close(outCon), add = TRUE)
        for (fastaIn in c(transcriptomeFasta, genomeFasta)) {
            inCon <- if (endsWith(fastaIn, ".gz")) {
                gzcon(file(fastaIn, open = "rb"))
            } else {
                file(fastaIn, open = "r")
            }
            repeat {
                chunk <- readLines(con = inCon, n = 10000L, warn = FALSE)
                if (length(chunk) == 0L) {
                    break
                }
                writeLines(text = chunk, con = outCon)
            }
            close(inCon)
        }
    }
    out <- c(gentrome = gentromeFile, decoys = decoyFile)
    assert(allAreFiles(out))
    alertSuccess("Salmon decoy files written.")
    invisible(out)
}
