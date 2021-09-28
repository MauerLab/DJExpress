.onLoad = function (libname, pkgname) {
  in.file = system.file("extdata", "junct.quant", package = "DJExpress")
  assign("in.file", in.file, envir = topenv())

  gtf0 = system.file("extdata", "chr1.gtf.gz", package = "DJExpress")
  assign('gtf0', gtf0, envir = topenv())
}
