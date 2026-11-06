# Build the local backup gene annotation used by process_uploaded_files.R when BioMart is unavailable.
# Run from the RNAview folder on a machine that can reach Ensembl:
#   Rscript gene_annotation_prep.R                                       # release 116 (tries e116, then the jun2026 archive)
#   Rscript gene_annotation_prep.R https://oct2024.archive.ensembl.org   # an older release from its Ensembl archive site
# Writes db/GeneAnnotation_human.rds, db/GeneAnnotation_mouse.rds and db/GeneAnnotation_rat.rds,
# with the same columns the upload page queries from BioMart.

#query the BioMart martservice of one Ensembl site directly (no Ensembl REST lookups as in biomaRt::useEnsembl),
#filter values are sent in chunks of 500 like biomaRt::getBM; returns a data.frame
martQuery <- function(host, dataset, attributes, filters=NULL, values=NULL) {
  url <- paste0(host, "/biomart/martservice")
  queryXML <- function(one) {
    flt <- if (is.null(one)) "" else sprintf('<Filter name="%s" value="%s"/>', filters, paste(one, collapse=","))
    sprintf('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" completionStamp="1" datasetConfigVersion="0.6"><Dataset name="%s" interface="default">%s%s</Dataset></Query>',
            dataset, flt, paste(sprintf('<Attribute name="%s"/>', attributes), collapse=""))
  }
  readResp <- function(resp) {
    httr::stop_for_status(resp)
    lines <- strsplit(httr::content(resp, as="text", encoding="UTF-8"), "\n", fixed=TRUE)[[1]]
    if (length(lines)==0 || lines[length(lines)]!="[success]") stop("incomplete BioMart response: ", substr(paste(lines, collapse=" "), 1, 200))
    read.delim(text=paste(lines[-length(lines)], collapse="\n"), header=TRUE, sep="\t", quote="", comment.char="",
               colClasses="character", na.strings="", check.names=FALSE, col.names=attributes)
  }
  #some archive sites lose a POSTed query behind a redirect (HTTP 405, or the martview page comes back),
  #so the query is sent again as GET ?query=..., with fewer filter values per request to keep the URL short
  getQuery <- function(one) {
    sub <- if (is.null(one)) list(NULL) else split(one, ceiling(seq_along(one)/100))
    do.call(rbind, lapply(sub, function(s) readResp(httr::GET(url, query=list(query=queryXML(s)), httr::timeout(600)))))
  }
  chunks <- if (is.null(filters)) list(NULL) else split(unique(values), ceiling(seq_along(unique(values))/500))
  usePost <- TRUE
  res <- lapply(chunks, function(one) {
    if (usePost) {
      out <- tryCatch(readResp(httr::POST(url, body=list(query=queryXML(one)), encode="form", httr::timeout(600))),
                      error=function(e) {cat("BioMart POST query failed on", host, ":", conditionMessage(e), "- trying GET\n"); NULL})
      if (!is.null(out)) return(out)
      usePost <<- FALSE
    }
    getQuery(one)
  })
  res <- do.call(rbind, res)
  if ("entrezgene_id" %in% names(res)) res$entrezgene_id <- as.integer(res$entrezgene_id)
  return(res)
}

args <- commandArgs(trailingOnly=TRUE)
datasets <- c(human="hsapiens_gene_ensembl", mouse="mmusculus_gene_ensembl", rat="rnorvegicus_gene_ensembl")
E_attributes <- c("ensembl_gene_id", "external_gene_name", "gene_biotype", "entrezgene_id", "description", "chromosome_name")

# export all three species from one Ensembl site, so they share a release
exportAll <- function(host) {
  # the release is only recorded with the files: when the registry does not report it, export anyway
  release <- tryCatch({
    txt <- httr::content(httr::GET(paste0(host, "/biomart/martservice?type=registry"), httr::timeout(60)), as="text", encoding="UTF-8")
    regmatches(txt, regexpr("(?<=Ensembl Genes )[0-9]+", txt, perl=TRUE))
  }, error=function(e) character(0))
  if (length(release)==0) {
    release <- paste0("unknown (exported ", Sys.Date(), ")")
    cat("Warning: cannot read the Ensembl release from", host, "\n")
  }
  # stop before downloading when the local backup is from a newer release than this site serves
  newRelease <- suppressWarnings(as.numeric(release))
  if (is.finite(localRelease) && is.na(newRelease)) stop("cannot confirm the release is not older than the local backup (Ensembl release ", localRelease, ")")
  if (is.finite(localRelease) && newRelease < localRelease) {
    stop(structure(class=c("olderRelease", "error", "condition"), list(call=NULL, message=paste0(
      "the local backup annotation (Ensembl release ", localRelease, ") is newer than the release on ", host, " (", release,
      "), so the files in db/ were kept. To replace them anyway, delete db/GeneAnnotation_*.rds and run this script again."))))
  }
  res <- list()
  for (sp in names(datasets)) {
    cat("Exporting", sp, "genes from Ensembl release", release, "on", host, "...\n")
    res[[sp]] <- martQuery(host, datasets[sp], E_attributes)
    attr(res[[sp]], "ensembl_release") <- release
  }
  return(res)
}

# newest Ensembl release among the existing local backup files (-Inf when there are none, or none records a numeric release)
localRelease <- suppressWarnings(max(sapply(names(datasets), function(sp) {
  strF <- file.path("db", paste0("GeneAnnotation_", sp, ".rds"))
  if (!file.exists(strF)) return(NA)
  as.numeric(attr(readRDS(strF), "ensembl_release"))
}), na.rm=TRUE))
if (is.finite(localRelease)) cat("Local backup annotation is Ensembl release", localRelease, "\n")

# www.ensembl.org moved to the new Ensembl platform (Aug 2026) without BioMart; BioMart is frozen at release 116 on the archive sites (e116, also served as the jun2026 archive); the useast/asia mirrors are retired
hosts <- if (length(args)>0) args[1] else c("https://e116.ensembl.org", "https://jun2026.archive.ensembl.org")
annots <- NULL
for (host in hosts) {
  annots <- tryCatch(exportAll(host), olderRelease=function(e) e, error=function(e) {
    cat("BioMart export failed on", host, ":", conditionMessage(e), "\n")
    NULL
  })
  if (inherits(annots, "olderRelease")) stop(conditionMessage(annots), call.=FALSE)
  if (!is.null(annots)) break
}
if (is.null(annots)) stop("BioMart is unavailable on all Ensembl sites, please try again later.")

dir.create("db", showWarnings=FALSE)
for (sp in names(annots)) {
  strF <- file.path("db", paste0("GeneAnnotation_", sp, ".rds"))
  saveRDS(annots[[sp]], strF, compress="xz")
  cat(nrow(annots[[sp]]), "rows,", length(unique(annots[[sp]]$ensembl_gene_id)), "genes saved to", strF, "\n")
}
