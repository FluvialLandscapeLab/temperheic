#' @export
thSeries = function(sg, xVals, tVals, specificUnits, ayeScale = 1, beeScale = 1) {
  if(!identical(attr(sg, "specificUnits"), specificUnits)) stop ("The units of the signal are not the same as those specified in specificUnits.")

  # make variables from required attributes.
  assign("aq", attr(sg, "aquifer"))
  assign("bd", attr(sg, "boundary"))
  assign("hy", attr(sg, "hydro"))

  # variable in eqn 8, Luce et al 2013
  a = (1 / (2 * hy$diffusivity_effective)) * (sqrt(((sqrt((hy$advectiveThermVel^4) + (4 * 2 * pi * bd$frequency * hy$diffusivity_effective)^2) + hy$advectiveThermVel^2)/2)) - hy$advectiveThermVel)
  b = (1 / (2 * hy$diffusivity_effective)) * sqrt(((sqrt((hy$advectiveThermVel^4) + (4 * 2 * pi * bd$frequency * hy$diffusivity_effective)^2) - hy$advectiveThermVel^2)/2))

  a = a * ayeScale
  b = b * beeScale

  # if (ayeScale > 0 & beeScale > 0)
  #   seriesList =
  #   lapply(
  #     xVals,
  #     function(x) bd$mean + bd$amplitude * exp(-((ayeScale * (x/1000) * a)) * x) * cos(((2 * pi / bd$period) * tVals - ((beeScale * (x/1000) * b)) * x) - (2 * pi * bd$phase / bd$period))
  #   )
  # else if (ayeScale > 0)
  #   seriesList =
  #   lapply(
  #     xVals,
  #     function(x) bd$mean + bd$amplitude * exp(-((ayeScale * (x/1000) * a)) * x) * cos((2 * pi / bd$period) * tVals - (b * x) - (2 * pi * bd$phase / bd$period))
  #   )
  # else if (beeScale > 0)
  #   seriesList =
  #   lapply(
  #     xVals,
  #     function(x) bd$mean + bd$amplitude * exp(-a * x) * cos(((2 * pi / bd$period) * tVals - ((beeScale * (x/1000) * b)) * x) - (2 * pi * bd$phase / bd$period))
  #   )
  # else
  seriesList =
    lapply(
      xVals,
      function(x) bd$mean + bd$amplitude * exp(-a * x) * cos((2 * pi / bd$period) * tVals - (b * x) - (2 * pi * bd$phase / bd$period))
    )
  names(seriesList) = paste0("x", xVals)

  # if (ayeScale > 0)
  #   amplitude =
  #   sapply(
  #     xVals,
  #     function(x)  bd$amplitude * exp(-((ayeScale * (x/1000) * a) *x)
  #     )
  #   )
  # else

  amplitude =
    sapply(
      xVals,
      function(x)  bd$amplitude * exp(-a *x)
    )
  names(amplitude) = names(seriesList)

  #   if (beeScale > 0)
  #   phaseRadians =
  #   sapply(
  #     xVals,
  #     function(x) (beeScale * (x/1000) * b) * x + (2 * pi * bd$phase / bd$period)
  #   )
  # else

  phaseRadians =
    sapply(
      xVals,
      function(x) (b * x + (2 * pi * bd$phase / bd$period))
    )
  names(phaseRadians) = names(seriesList)

  factorialDists = expand.grid(from = names(seriesList), to = names(seriesList))

  ampRatio = amplitude[factorialDists$to] / amplitude[factorialDists$from]
  deltaPhaseRadians = phaseRadians[factorialDists$to] - phaseRadians[factorialDists$from]
  eta = -log(ampRatio) / deltaPhaseRadians
  names(xVals) = names(seriesList)
  deltaXvals = (xVals[factorialDists$to] - xVals[factorialDists$from])
  diffusivity_effectiveEmp = (eta * (2*pi / bd$period) * deltaXvals^2) / (log(ampRatio)^2 + deltaPhaseRadians^2)
  advectiveThermVelEmp = (((2*pi / bd$period) * deltaXvals) / sqrt((log(ampRatio)^2 + deltaPhaseRadians^2))) * ((1 - eta^2) / sqrt(1 + eta^2))

  derivedVals = array(
    data = c(ampRatio, deltaPhaseRadians, eta),
    dim = c(length(xVals), length(xVals), 3),
    dimnames = list(from = names(seriesList), to = names(seriesList), value = c("ampRatio", "deltaPhaseRadians", "eta"))
  )

  empiricalVals = array(
    data = c(diffusivity_effectiveEmp, advectiveThermVelEmp),
    dim = c(length(xVals), length(xVals), 2),
    dimnames = list(from = names(seriesList), to = names(seriesList), value = c("diffusivity_effectiveEmp", "advectiveThermVelEmp"))
  )

  DF = structure(
    do.call(data.frame, args = c(list(time = tVals), seriesList)),
    signal = sg,
    phaseRadians = phaseRadians,
    amplitude = amplitude,
    derivedVals = derivedVals,
    empiricalVals = empiricalVals,
    class = c("thSeries", "data.frame")
  )

  return(DF)
}

# #' @export
# as.xts.thSeries = function(x, POSIXctRange, xRange = c(1, nrow(x)), trim = F, ...) {
#   if(trim) {
#     x = x[xRange[1]:xRange[2],]
#   }
#   numericRange = as.numeric(POSIXctRange)
#   time = x$time[xRange]
#   time2Date = lm(POSIXctRange ~ time)
#   POSIXctIndex = as.POSIXct(round(predict(time2Date, newdata = x)), origin = "1970-01-01 00:00.00 UTC")
#   class(x) = "data.frame"
#   return(xts::as.xts(x, order.by = POSIXctIndex))
# }
