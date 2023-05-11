library(tidyverse)



Creat_curve <- function(Temp, dH, dCp, Tm) {
  Cp <- (((dH+(dCp*(Temp-Tm)))*((dH/Temp)
                                   -((dCp*Tm)/Temp)+dCp)
                    *(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                          /(0.00198588*Temp))/(0.00198588*Temp))
                    /((exp(-((dH*(1-(Temp/Tm)))
                             +(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))/(0.00198588*Temp))+1)^2))
                   +((dCp*(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                               /(0.00198588*Temp))))
                     /(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                           /(0.00198588*Temp))+1)))
  return (Cp)
}



