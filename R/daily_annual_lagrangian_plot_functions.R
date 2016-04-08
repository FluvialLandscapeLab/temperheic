daily_temp_plotter = function(x, distance, lineColor, lineType = 1:1000){
  disp.L = dispersivity
  mean.temp = 22.02152# based on mean 2012 stream temp from fit C
  stream_amp = 3.063285  # based on amplitude 2012 stream temp from fit C
  stream_phase = 3.340407 # based on phase 2012 stream temp from fit C
  w.annual = (2 * pi) / (86400)
  year.seconds = seq(10 * 3600, (70 * 3600 - 3600), 3600)
  times = 1:48 #seq(ymd_hms("2012-7-1 10:00:00"), ymd_hms("2012-7-4 09:00:00"), 3600)
  D = (x$disp.L * x$v.x) ## mechanical dispersion (m2/s see Zheng 1999)
  ke = (kt.bulk / c.vol)  + (disp.L) * (((p.w * cw ) / c.vol) * q)
  # exponent in eqn 8, Luce et al 2013
  a = (1/(2*ke))*(sqrt(((sqrt((vt^4) + (4*w.annual*ke)^2) + vt^2)/2))-vt)
  # variable in eqn 8, Luce et al 2013
  b = (1/(2*ke))* sqrt(((sqrt((vt^4) + (4*w.annual*ke)^2) - vt^2)/2))
  T.mu = mean.temp
  z = distance #seq(0,5,.1)
  #    T.LG.1.D.1 = T.mu + stream_amp*exp(-a[1]*z[1])*cos((w.annual * (year.seconds) - b[1]*z[1]) - stream_phase)
  #    T.LG.1.D.2 = T.mu + stream_amp*exp(-a[2]*z[2])*cos((w.annual * (year.seconds) - b[2]*z[2]) - stream_phase)
  #    T.LG.1.D.3 = T.mu + stream_amp*exp(-a[3]*z[3])*cos((w.annual * (year.seconds) - b[3]*z[3]) - stream_phase)
  #
  #    T.LG.1.D.1.TS = zoo(T.LG.1.D.1, order.by = times)
  #    T.LG.1.D.2.TS = zoo(T.LG.1.D.2, order.by = times)
  #    T.LG.1.D.3.TS = zoo(T.LG.1.D.3, order.by = times)
  
  lineColor = rep(lineColor, length(z)/length(lineColor))
  
  for (i in 1:length(z)) {
    T.LG.1.D = T.mu + stream_amp*exp(-a[i]*z[i])*cos((w.annual * (year.seconds) - b[i]*z[i]) - stream_phase)
    T.LG.1.D.TS = zoo(T.LG.1.D, order.by = times)
    lines(T.LG.1.D.TS, lwd = 1, col = lineColor[i], lty = lineType[i])
  }
}


#############  function to create lines for vectors of dispersivity
daily_temp_plotter = function(dispersivity, distance, lineColor, lineType = 1:1000){
    disp.L = dispersivity
    mean.temp = 22.02152# based on mean 2012 stream temp from fit C
    stream_amp = 3.063285  # based on amplitude 2012 stream temp from fit C
    stream_phase = 3.340407 # based on phase 2012 stream temp from fit C
    w.annual = (2 * pi) / (86400)
    year.seconds = seq(10 * 3600, (70 * 3600 - 3600), 3600)
    times = 1:48 #seq(ymd_hms("2012-7-1 10:00:00"), ymd_hms("2012-7-4 09:00:00"), 3600)
    D = (disp.L * v.x) ## mechanical dispersion (m2/s see Zheng 1999)
    ke = (kt.bulk / c.vol)  + (disp.L) * (((p.w * cw ) / c.vol) * q)
    # exponent in eqn 8, Luce et al 2013
    a = (1/(2*ke))*(sqrt(((sqrt((vt^4) + (4*w.annual*ke)^2) + vt^2)/2))-vt)
    # variable in eqn 8, Luce et al 2013
    b = (1/(2*ke))* sqrt(((sqrt((vt^4) + (4*w.annual*ke)^2) - vt^2)/2))
    T.mu = mean.temp
    z = distance #seq(0,5,.1)
#    T.LG.1.D.1 = T.mu + stream_amp*exp(-a[1]*z[1])*cos((w.annual * (year.seconds) - b[1]*z[1]) - stream_phase)
#    T.LG.1.D.2 = T.mu + stream_amp*exp(-a[2]*z[2])*cos((w.annual * (year.seconds) - b[2]*z[2]) - stream_phase)
#    T.LG.1.D.3 = T.mu + stream_amp*exp(-a[3]*z[3])*cos((w.annual * (year.seconds) - b[3]*z[3]) - stream_phase)
#
#    T.LG.1.D.1.TS = zoo(T.LG.1.D.1, order.by = times)
#    T.LG.1.D.2.TS = zoo(T.LG.1.D.2, order.by = times)
#    T.LG.1.D.3.TS = zoo(T.LG.1.D.3, order.by = times)
    
    lineColor = rep(lineColor, length(z)/length(lineColor))
    
    for (i in 1:length(z)) {
      T.LG.1.D = T.mu + stream_amp*exp(-a[i]*z[i])*cos((w.annual * (year.seconds) - b[i]*z[i]) - stream_phase)
      T.LG.1.D.TS = zoo(T.LG.1.D, order.by = times)
      lines(T.LG.1.D.TS, lwd = 1, col = lineColor[i], lty = lineType[i])
    }
    
    
    # vd = sqrt(ke*w.annual*2)
    # hourLag = (c(1,2,3)/vt)/3600 #(c(1,2,3)/vd)/3600
    # hourLagFloor = hourLag%/%1
    # 
    # timeOffset1 = hourLagFloor[1]:(24+hourLagFloor[1])
    # timeOffset2 = hourLagFloor[2]:(24+hourLagFloor[2])
    # timeOffset3 = hourLagFloor[3]:(24+hourLagFloor[3])
    # #print(timeOffset1); print(timeOffset2); print(timeOffset3)
    ########## now print lines for daily temps
    #####
    #[hourLagFloor[1]:(24+hourLagFloor[1])]
    # lines(T.LG.1.D.1.TS, lwd = 1, col = lineColor, lty = 1)
    # lines(T.LG.1.D.2.TS, lwd = 1, col = lineColor, lty = 2)
    # lines(T.LG.1.D.3.TS, lwd = 1, col = lineColor, lty = 3)

    #text(, , bquote(beta == .(disp.L)), pos = 4, cex = .7, col = "black")
}

daily_input_temp_plotter = function(dispersivity, distance, lineColor){
  disp.L = dispersivity
  mean.temp = 22.02152# based on mean 2012 stream temp from fit C
  stream_amp = 3.063285  # based on amplitude 2012 stream temp from fit C
  stream_phase = 3.340407 # based on phase 2012 stream temp from fit C
  w.annual = (2 * pi) / (86400)
  year.seconds = seq(10 * 3600, (34 * 3600 - 3600), 3600)
  D = (disp.L * v.x) ## mechanical dispersion (m2/s see Zheng 1999)
  ke = (kt.bulk / c.vol)  + (disp.L) * (((p.w * cw ) / c.vol) * q)
  # exponent in eqn 8, Luce et al 2013
  a = (1/(2*ke))*(sqrt(((sqrt((vt^4) + (4*w.annual*ke)^2) + vt^2)/2))-vt)
  # variable in eqn 8, Luce et al 2013
  b = (1/(2*ke))* sqrt(((sqrt((vt^4) + (4*w.annual*ke)^2) - vt^2)/2))
  T.mu = mean.temp
  z = distance #seq(0,5,.1)
  T.LG.1.D.1 = T.mu + stream_amp*exp(-a[1]*z)*cos((w.annual * (year.seconds) - b[1]*z) - stream_phase)
  T.LG.1.D.2 = T.mu + stream_amp*exp(-a[2]*z)*cos((w.annual * (year.seconds) - b[2]*z) - stream_phase)
  T.LG.1.D.3 = T.mu + stream_amp*exp(-a[3]*z)*cos((w.annual * (year.seconds) - b[3]*z) - stream_phase)

  ########## now print lines for daily temps
  #####
  lines(T.LG.1.D.1, lwd = 1, col = lineColor, lty = 1)
  
  #text(, , bquote(beta == .(disp.L)), pos = 4, cex = .7, col = "black")
}


annual_temp_plotter = function(dispersivity, distance, lineColor){
    disp.L = dispersivity
    mean.temp = 11.34786# based on mean 2012 stream temp from fit C
    stream_amp = 10.38095  # based on amplitude 2012 stream temp from fit C
    stream_phase = 0.37 # based on phase 2012 stream temp from fit C
    w.annual = (2 * pi) / (86400 * 366) # frequency
    year.seconds = seq(0, (366 * 86400 - 86400), 86400)
    D = (disp.L * v.x) ## mechanical dispersion (m2/s see Zheng 1999)
    ke = (kt.bulk / c.vol)  + (disp.L) * (((p.w * cw ) / c.vol) * q)
    # exponent in eqn 8, Luce et al 2013
    a = (1/(2*ke))*(sqrt(((sqrt((vt^4) + (4*w.annual*ke)^2) + vt^2)/2))-vt)
    # variable in eqn 8, Luce et al 2013
    b = (1/(2*ke))* sqrt(((sqrt((vt^4) + (4*w.annual*ke)^2) - vt^2)/2))
    T.mu = mean.temp
    z = distance #seq(0,2000,5)
    T.LG.1.D.1 = T.mu + stream_amp*exp(-a[1]*z[1])*cos((w.annual * (year.seconds) - b[1]*z[1]) - stream_phase)
    T.LG.1.D.2 = T.mu + stream_amp*exp(-a[2]*z[2])*cos((w.annual * (year.seconds) - b[2]*z[2]) - stream_phase)
    T.LG.1.D.3 = T.mu + stream_amp*exp(-a[3]*z[3])*cos((w.annual * (year.seconds) - b[3]*z[3]) - stream_phase)
    
    ########## now print lines for annual temps
    #####

    lines(T.LG.1.D.1, lwd = 1, col = lineColor, lty = 1)
    lines(T.LG.1.D.2, lwd = 1, col = lineColor, lty = 2)
    lines(T.LG.1.D.3, lwd = 1, col = lineColor, lty = 3)
    #text(, , bquote(beta == .(disp.L)), pos = 4, cex = .7, col = "black")
}
