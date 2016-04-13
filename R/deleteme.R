

#' @export
as.xts.thSeries = function(x, POSIXct.origin) {
  theOrder = as.POSIXct(x$time, origin = POSIXct.origin)
  x = x[-match("time", names(x))]
  x = xts(x, order.by = theOrder)
  return(x)
}
