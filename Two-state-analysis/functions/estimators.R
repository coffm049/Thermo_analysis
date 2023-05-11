####################################################################################################################
########## Step 4. Graphical features ( raw dCp, dH, and Tm)   #####################################################
####################################################################################################################
e     <- 2.718281828459
R     <- 0.00198588
T.end <- 338 # Choose a temperature at which the transition has ended, but not two far away.

#' Estimate thermodynamic features from geometrical data of 2-state calorimetric process
#'
#' @param Temp vector of temperatures (K).
#' @param Cp vector of heat capacities (kcal/mol T).
#' @param T.begin temperature at which the transition begins.
#' @param T.end temperature at which the transition ends.
#' @return vals a the estimated dH, dCp, and Tm, and Vant hoff estimated enthalpy.
#' @examples
#' data(lysozyme)
#' geom_ests(lysozyme$Temp, lysozyme$Cp, 330, 390)
geom_ests <- function(Temp, Cp, T.begin, T.end) {
  # Choose index Temp closest to T.end
  T.end.index  <- which(abs(Temp-T.end)==min(abs(Temp-T.end)))
  # dCp from taking difference between baselines
  dCp.geom     <- mean(tail(Cp,3))-mean(head(Cp,3))
# sets the lower bound for integration based on user input for end of folded baseline
  int.l        <-  which(abs(Temp-T.begin)==min(abs(Temp-T.begin)))
  int.data <- cumsum(Cp[int.l:T.end.index]*diff(Temp[int.l:(T.end.index+1)]))
  signal.sigmoid <- int.data/tail(int.data,n=1) * dCp.geom
  trans.sig <- Cp[int.l:T.end.index]-signal.sigmoid
  int.min.sig <- cumsum(trans.sig*diff(Temp[int.l:(T.end.index+1)]))
  norm.int.min.sig <- int.min.sig/tail(int.min.sig,1)
  dH.geom <- tail(int.min.sig, n=1)
  dH.half.geom <- dH.geom/2
  Tm.geom.index  <- which(abs(dH.half.geom-int.min.sig)==min(abs(dH.half.geom-int.min.sig)))+int.l
  Tm.geom <- Temp[Tm.geom.index]

  # Vant hoff estimation
  dH.VH              <- 4*R*(Tm.geom^2)*mean((diff(norm.int.min.sig[(Tm.geom-10):(Tm.geom+10)])/diff(Temp[(Tm.geom-10):(Tm.geom+10)])))
  return(c(dH.geom, Tm.geom, dCp.geom, dH.VH))
}



#' Estimate thermodyanmic parameters of calorimetric two-state transition with two state nonlinear model
#' 
two.state <- function(df, guess_Tm = 290, guess_dCp = 3, guess_dH = 100) {
  fit.params <- nls(Cp~
                      (((dH+(dCp*(Temp-Tm)))*((dH/Temp)
                                              -((dCp*Tm)/Temp)+dCp)
                        *(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                              /(0.00198588*Temp))/(0.00198588*Temp))
                        /((exp(-((dH*(1-(Temp/Tm)))
                                 +(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))/(0.00198588*Temp))+1)^2))
                       +((dCp*(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                                   /(0.00198588*Temp))))
                         /(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                               /(0.00198588*Temp))+1))),
                    start=list(Tm=guess_Tm, dCp=guess_dCp, dH=guess_dH),
  data = df)
  return(fit.params)
}


