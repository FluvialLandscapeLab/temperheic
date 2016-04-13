#' @export
thSeries = function(sg, xVals, tVals, specificUnits) {
  if(!identical(attr(sg, "specificUnits"), specificUnits)) stop ("The units of the signal are not the same as those specified in specificUnits.")

  # make variables from required attributes.
  assign("aq", attr(sg, "aquifer"))
  assign("bd", attr(sg, "boundary"))
  assign("hy", attr(sg, "hydro"))

  # variable in eqn 8, Luce et al 2013
  a = (1/(2*hy$diffusivity_effective))*(sqrt(((sqrt((hy$advectiveThermVel^4) + (4*bd$frequency*hy$diffusivity_effective)^2) + hy$advectiveThermVel^2)/2))-hy$advectiveThermVel)
  b = (1/(2*hy$diffusivity_effective))* sqrt(((sqrt((hy$advectiveThermVel^4) + (4*bd$frequency*hy$diffusivity_effective)^2) - hy$advectiveThermVel^2)/2))

  seriesList =
    lapply(
      xVals,
      function(x) bd$mean + bd$amplitude*exp(-a*x)*cos(((2*pi/bd$period) * tVals - b*x) - (2*pi*bd$phase/bd$period))
    )
  names(seriesList) = paste0("x", xVals)
  return(
    structure(
      do.call(data.frame, args = c(list(time = tVals), seriesList)),
      class = c("thSeries", "data.frame")
    )
  )
}

#' @export
as.xts.thSeries = function(x, POSIXctRange, xRange = c(1, nrow(x)), trim = F, ...) {
  if(trim) {
    x = x[xRange[1]:xRange[2],]
  }
  numericRange = as.numeric(POSIXctRange)
  time = x$time[xRange]
  time2Date = lm(POSIXctRange ~ time)
  POSIXctIndex = as.POSIXct(round(predict(time2Date, newdata = x)), origin = "1970-01-01 00:00.00 UTC")
  class(x) = "data.frame"
  return(xts::as.xts(x, order.by = POSIXctIndex))
}
