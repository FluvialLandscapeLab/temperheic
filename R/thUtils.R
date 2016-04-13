
# Replaces general units (e.g., "E t-1 L-3 T-1") with specific units specified
# by the user (e.g., "kJ s-1 m-3 degC-1")
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
