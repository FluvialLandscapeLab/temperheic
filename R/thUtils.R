#' @export
as.xts.thSeries = function(x, POSIXct.origin) {
  theOrder = as.POSIXct(x$time, origin = POSIXct.origin)
  x = x[-match("time", names(x))]
  x = xts(x, order.by = theOrder)
  return(x)
}

#' @export
as.zoo.thSeries = function(x, POSIXct.origin) {
  theOrder = as.POSIXct(x$time, origin = POSIXct.origin)
  x = x[-match("time", names(x))]
  x = zoo(x, order.by = theOrder)
  return(x)
}

#' @export
htPlot = function(myThSeries, POSIXct.origin = "2014-01-01 00:00:00") {
  myXTS = as.xts.thSeries(myThSeries, POSIXct.origin)
  plot(myXTS$x0, main = "", ylab = "Temperature (degC)", type = 'n', auto.grid = F, minor.ticks = F)
  lines(myXTS$x0, col = "grey70", lwd = 3, lty=2)
  for(col in names(myXTS)[2:length(names(myXTS))]) {
    lines(myXTS[,col], lwd = 3)
  }
  seriesAttr = attributes(myThSeries)
  sg = attr(myThSeries, "signal")
  hy = attr(sg, "hydro")
  k = hy$hydCond
  d = hy$dispersivity
  text(as.numeric(min(index(myXTS))), max(myXTS$x0), pos = 4, paste0("k = ", round(k,6), "  B = ", round(d, 6)))
}


# Replaces general units (e.g., "E t-1 L-3 T-1") with specific units specified
# by the user (e.g., "kJ s-1 m-3 degC-1")
#' @export
.thSpecificUnits = function(generalUnits, specificUnits = thUnits()) {
  atomicUnits = unlist(lapply(generalUnits, strsplit, split = " "), recursive = F)

  #regular expression for an optional "-" and any number of digits at the end of
  #a string
  pattern = "[-]?[[:digit:]]+$"
  #strip off any match to pattern
  generalSymbols = lapply(atomicUnits, sub, pattern = pattern, replacement = "")
  #return any substring that matches to pattern
  generalExponents =
    lapply(
      atomicUnits,
      function(x) {
        locations = regexpr(pattern = pattern, x)
        locations[locations == -1] = 1000000L
        substring(x, locations)
      }
    )

  #return any generalSymbols that are not expected.  If any, throw error and
  #report to user
  unexpectedGenUnits = sapply(generalSymbols, function(x) {any(!(x %in% names(specificUnits)))})
  if(any(unexpectedGenUnits)) {
    cat("Unexpected general units found in: ", paste0("'", generalUnits[unexpectedGenUnits], "'", collapse = ", "), '\n  Expected units are: ', paste0(attributes(specificUnits)$longUnitName, "='", names(specificUnits), "'", collapse = ", "))
    stop("A general unit can be followed by positive or negative number (representing and exponent).\nThere can be no space between a unit and its exponent (e.g., 'M-3' not 'M -3').\nThere must be spaces between unit/exponent pairs (e.g. 'E L-3', not 'EL-3')")
  }

  #replace the generalSymbols with corresponding specific ones
  generalSymbols = lapply(generalSymbols, function(x) unlist(specificUnits)[x])
  #concatinate specific units with exponents and return vector specific units and exponents
  return(mapply(paste0, generalSymbols, generalExponents, collapse = " "))
}
