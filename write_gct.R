write_gct <- function(ds, ofile, precision=4, appenddim=TRUE, ver=3) {
  if (!methods::is(ds, "GCT")) {
    stop("ds must be a GCT object")
  }
  # make sure it's valid
  methods::validObject(ds)
  
  # extract the components
  m <- mat(ds)
  rdesc <- meta(ds)
  cdesc <- meta(ds, dimension="column")
  rid <- ids(ds)
  cid <- ids(ds, dimension="column")
  
  # append the dimensions of the data set, if desired
  if (appenddim) ofile <- append_dim(ofile, m, extension="gct")
  
  precision <- floor(precision)
  cat("Saving file to ", ofile, "\n")
  nr <- nrow(m)
  nc <- ncol(m)
  cat(sprintf("Dimensions of matrix: [%dx%d]\n", nr, nc))
  cat(sprintf("Setting precision to %d\n", precision))
  # open file and write   
  if (ver == 3) {
    # remove the 'id' columns
    cdesc$id <- NULL
    rdesc$id <- NULL
    # get the counts of meta data fields
    nrdesc <- ncol(rdesc)
    ncdesc <- ncol(cdesc)
    colkeys <- names(cdesc)
    # append header
    cat(sprintf("#1.%d\n%d\t%d\t%d\t%d", ver, nr, nc, nrdesc, ncdesc),
        file=ofile,sep='\n')      
    # line 3: sample row desc keys and sample names
    cat(paste(c("id", names(rdesc), cid), collapse="\t"),
        file=ofile, sep="\n", append=TRUE)
    # line 4 + ncdesc: sample desc
    filler <- 'na'
    if (ncdesc > 0) {
      for (ii in seq_len(ncdesc)) {
        if (is.numeric(cdesc[, ii])) {
          cat(paste(c(colkeys[ii], rep(filler, nrdesc),
                      round(cdesc[, ii], precision)),
                    collapse="\t"),
              file=ofile, sep="\n", append=TRUE)  
        } else {
          cat(paste(c(colkeys[ii], rep(filler, nrdesc),
                      cdesc[, ii]),
                    collapse="\t"),
              file=ofile, sep="\n", append=TRUE)
        }
      }
    }
    for (ii in seq_len(nr)) {    
      # print rows
      cat(paste(c(rid[ii],
                  rdesc[ii, ],
                  round(m[ii, ], precision)), collapse="\t"),
          sep="\n", file=ofile, append=TRUE)
    }
  } else {
    # assume ver 1.2 and below, ignore descriptors
    # append header
    cat(sprintf("#1.%d\n%d\t%d", ver, nr, nc),
        file=ofile, sep="\n")      
    # line 3: sample row desc keys and sample names
    cat(paste(c("id", "Description", cid), collapse="\t"),
        file=ofile, sep="\n", append=TRUE)
    for (ii in seq_len(nr)) {    
      # print rows
      cat(paste(c(rid[ii],
                  rdesc[ii, 2],
                  round(m[ii, ], precision)), collapse="\t"),
          sep="\n", file=ofile, append=TRUE)
    }
  }
  cat("Saved.\n")  
}