setwd("/home/christian/Desktop")  # defines where the files are.
# NOTE: this file path will not work for you this is just an example. So you need to insert your own file path in between the quotations above.
# Before reading data into R, make sure you saved it as a csv file.
# -----------------------------------------------------------------------------------------------------------#
########## step 1. Reading in and plotting raw Data ##########################################################
#This requires all Cp values from multiple scans (if multiple scans are being analyzed) to be in the same spreadsheet.
##############################################################################################################
df.dsc  <- read.csv(file.choose(), header = TRUE, sep=","); # defines the file beng used
#The value after $ is the name of the column you wish to read. So the cp of scan 1 was called "X1cp" in the spreadsheet"
#Same thing for the next set of cp columns. If there aren't any more values you wish to read, you can skip this step
# either by excepting the error values or commenting out this section of code.
if(sum(df.dsc[[1]]>270)>1) {df.dsc[[1]]} else {df.dsc[[1]] <-df.dsc[[1]]+273.15}

# Next, plot a raw data graph to make sure it looks as expected, pch=point size and type
par(mfrow=c(ceiling(sqrt(length(df.dsc)-1)), (length(df.dsc)-1)/ceiling(sqrt(length(df.dsc)-1))))
lapply(c(2:length(df.dsc)), function(i) {
  plot(x=df.dsc[[1]], y=df.dsc[[i]],
       main="DSC", xlab="Temp(K)", ylab="dCp (kcal/mol/K)", pch=".")
})
# At this point you should see a plot(s) of your raw data in the plots window on the right.
##############################################################################################################
########## step 2. Trimming off the initial low temperature fluctuation and plotting    ######################
##############################################################################################################
# Select a portion of the plot that is flat and will serve as your folded baseline. Enter the bounds for the temperature below.
t.lower       <- 305    ; trim.Lbound  <- which(abs(df.dsc[[1]]-t.lower)==min(abs(df.dsc[[1]]-t.lower)))    # the values entered do not need to match the values exactly. You can just right '310' for instance
t.upper       <- 310
trim.upper    <- 385    ; trim.Ubound  <- which(abs(df.dsc[[1]]-trim.upper)==min(abs(df.dsc[[1]]-trim.upper)))

# Plot the data with the fat trimmed off
par(mfrow=c(ceiling(sqrt(length(df.dsc)-1)), (length(df.dsc)-1)/ceiling(sqrt(length(df.dsc)-1))))
trim.data <- lapply(1:length(df.dsc), function (i) {
  plot(x=df.dsc[[1]][trim.Lbound:trim.Ubound], 
       y=df.dsc[[i]][trim.Lbound:trim.Ubound], 
       xlab="Temp(C)",ylab="Cp(kcal/mol/K)", main="Trim DSC")
  return(df.dsc[[i]][trim.Lbound:trim.Ubound])
}) ; trim.Lbound  <- which(abs(trim.data[[1]]-t.lower)==min(abs(trim.data[[1]]-t.lower))) ;trim.Ubound  <- which(abs(trim.data[[1]]-trim.upper)==min(abs(trim.data[[1]]-trim.upper)))
#############################################################################################################
########## step 3. subtracting off the baseline #############################################################
#############################################################################################################
trans.end          <- 337    # Choose a temperature at which the temperature is for sure over.
trans.end2         <- 342    # choose a temperature at which the unfolded baseline is still flat

par(mfrow=c((ceiling(sqrt(length(trim.data)-1))),2))
data.wo.base       <- lapply(2:length(trim.data), function(i) {
  t.end.1 <- which(abs(trim.data[[1]]-trans.end)==min(abs(trim.data[[1]]-trans.end)))
  t.end.2 <- which(abs(trim.data[[1]]-trans.end2)==min(abs(trim.data[[1]]-trans.end2))) 
  fbl.Ubound         <- which(abs(trim.data[[1]]-t.upper)==min(abs(trim.data[[1]]-t.upper)))
  fbl              <-  lm(trim.data[[i]][1:fbl.Ubound] ~ trim.data[[1]][1:fbl.Ubound])
  ubl              <-  lm(trim.data[[i]][t.end.1:t.end.2] ~ trim.data[[1]][t.end.1:t.end.2])
  plot(x=trim.data[[1]], y=trim.data[[i]],
       xlab="dCp(kcal/mol/K)",ylab="Temp(K)", main="DSC flat baselines", pch=".")
  lines(x=trim.data[[1]], y=(ubl[["coefficients"]][[1]]+ubl[["coefficients"]][[2]]*trim.data[[1]]), col='red')
  lines(x=trim.data[[1]], y=(fbl[["coefficients"]][[1]]+fbl[["coefficients"]][[2]]*trim.data[[1]]), col='red')
  plot(x=trim.data[[1]], y=(trim.data[[i]]-(fbl[["coefficients"]][[1]]+fbl[["coefficients"]][[2]]*trim.data[[1]])),
       xlab="dCp(kcal/mol/K)",ylab="Temp(K)", main="DSC wo baseline")
  return(trim.data[[i]]-(fbl[["coefficients"]][[1]]+fbl[["coefficients"]][[2]]*trim.data[[1]]))
})


########################END curvy baseline fit ########################################################
cp.total    <- lapply(c(1:length(data.wo.base[[1]])), function(i){
  return(sapply(c(1:length(data.wo.base)), function(j) {
    data.wo.base[[j]][i]
  }))
})
cp.ave <- sapply(c(1:length(cp.total)), function(i) {
  mean(cp.total[[i]])
})
cp.std <- sapply(c(1:length(cp.total)), function(i) {
  sd(cp.total[[i]])
})

ffbl.lm <- lm(cp.ave[trim.Lbound:trim.Ubound] ~ trim.data[[1]][trim.Lbound:trim.Ubound])                                      # line model for the baseline after baseline subtraction
ffbl <- (summary(ffbl.lm)$coefficients[1,1])+((summary(ffbl.lm)$coefficients[2,1])*trim.data[[1]])    # definition of said baseline

par(mfrow=c(1,1))
plot(trim.data[[1]]-273.15, data.wo.base[[1]],pch=".", xlab="Temp (C)" , ylab="Cp (kcal/mol C)",
     cex.axis=2, lwd=3, col='grey20')
lines(trim.data[[1]]-273.15, data.wo.base[[2]],pch=".", xlab="Temp (C)" , ylab="Cp (kcal/mol C)",
     cex.axis=2, lwd=3, col='grey55')  
lines(trim.data[[1]]-273.15, data.wo.base[[3]], col="blue")

# plots the new baseline. It should be flat.
lapply(c(1:length(cp.wo.baseline)), function(i) {
  lines(trim.data[[1]], cp.wo.baseline[[i]], col=i)
})
lines(trim.data[[1]], ffbl, lwd=4, col='red')
library("ggplot2")
df.plot <- data.frame(c(trim.data[[1]]))
df.plot$cp.ave <- cp.ave
colnames(df.plot)[1] <- "Temp"
ggplot(df.plot, aes(Temp, cp.ave)) + 
  geom_point()

####################################################################################################################
########## Step 4. Graphical features ( raw dCp, dH, and Tm)   #####################################################
####################################################################################################################
e     <- 2.718281828459
R     <- 0.00198588
T.end <- 338 # Choose a temperature at which the transition has ended, but not two far away.
par(mfrow=c(2,1), mar=c(4,7,2,1))
Graphical.features <- function(experiment) {
  T.end.index  <- which(abs(trim.data[[1]]-T.end)==min(abs(trim.data[[1]]-T.end)))                    # pick a temperature at which the transition has ended
  dCp.geom     <- mean(tail(experiment,40))-mean(head(experiment,40))                                 # dCp from taking difference between baselines
  int.l        <-  which(abs(trim.data[[1]]-t.upper)==min(abs(trim.data[[1]]-t.upper)))               # sets the lower bound for integration based on user input for end of folded baseline
  int.data     <- cumsum(experiment[int.l:T.end.index]*diff(trim.data[[1]][int.l:(T.end.index+1)]))
  plot (trim.data[[1]][int.l:T.end.index],int.data, ylab="Enthalpy (kcal/mol)", xlab= "Temp (K)")                                                  # plots the integrated change in heat capacity
  norm.int.data    <- int.data/tail(int.data,n=1)                                                     # normalizes the integral to the dH of transition
  plot (trim.data[[1]][int.l:T.end.index],norm.int.data, pch=".", ylab= "Normalized Enthalpy (a.u.)", xlab="Temp (K)")                                    # Plots the normalized integral
  signal.sigmoid   <- norm.int.data*dCp.geom                                                          # multiplies the normalized signal by the difference in heat capacity to create spline baseline
  plot(trim.data[[1]], experiment, pch=".", ylab="Heat capacity (kcal mol^-1 K^-1)", xlab="Temp (K)")                                                     
  lines(trim.data[[1]][int.l:T.end.index],signal.sigmoid)                                             # plot to see if spline does a good job for the baseline
  trans.sig          <- cp.ave[int.l:T.end.index]-signal.sigmoid                                      # subtracts off the spline baseline
  int.min.sig        <- cumsum(trans.sig*diff(trim.data[[1]][int.l:(T.end.index+1)]))                 # Integrates transition without spline baseline
  norm.int.min.sig   <- int.min.sig/tail(int.min.sig,1)                                               # normalized integral (i.e concentration of unfolded species)
  plot(trim.data[[1]][int.l:T.end.index], norm.int.min.sig, xlab = "Temp(K)", ylab="Extent of reaction")
  dH.geom            <- tail(int.min.sig, n=1)                                                        # determines dH of transion from total integral
  dH.half.geom       <- dH.geom/2                                                                     # determines half dH for Tm guess
  Tm.geom.index      <- which(abs(dH.half.geom-int.min.sig)==min(abs(dH.half.geom-int.min.sig)))+int.l
  Tm.geom            <- trim.data[[1]][Tm.geom.index]
  H.geom     <- dH.geom+(dCp.geom*(trim.data[[1]]-Tm.geom))    # define enthalpy using geom fit values
  G.geom     <- dH.geom*(1-(trim.data[[1]]/Tm.geom))+dCp.geom*(trim.data[[1]]-Tm.geom-trim.data[[1]]*log(trim.data[[1]]/Tm.geom))    # define free energy term
  K.geom     <- e^(-G.geom/(R*trim.data[[1]]))                                                        # defines  equilibrium
  dG37.geom.i<- which(abs(trim.data[[1]]-(273.15+37))==min(abs(trim.data[[1]]-(273.15+37))))
  dG37.geom  <- G.geom[dG37.geom.i]
  plot(trim.data[[1]]-273.15, experiment, pch=".",xlim= c(t.lower-273.15, T.end-271), ylim=c(-0.1, max(experiment)+1), 
       xlab=expression(paste("Temperature ( ",degree,"C)")) , 
       ylab= expression(paste("Cp (kcal/mol ",degree,"C)")), cex.axis=2, cex.lab=2)
  
  lines(trim.data[[1]][int.l:T.end.index]-273.15, signal.sigmoid)
  polygon(c(trim.data[[1]][int.l:T.end.index]-273.15, rev(trim.data[[1]][int.l:T.end.index])-273.15), c(experiment[int.l:T.end.index], rev(signal.sigmoid)), col=rgb(0.5,0,0.1,0.25), border=NA)
  text(x=Tm.geom-280, y=10, labels=expression(paste(Delta, "H")), cex=2, col='black')
  segments(x0=Tm.geom-280, y0=9.5, x1=Tm.geom-276.15, y1=4)
  segments(x0=Tm.geom-263, y0=0, x1=Tm.geom-263, y1=dCp.geom, lty=3)
  text(x=Tm.geom-267, y=dCp.geom/2, labels=expression(paste(Delta, "Cp")), cex=2, col='black')
  segments(x0=Tm.geom-273.15, y0=dCp.geom/2, x1=Tm.geom-273.15, y1=experiment[Tm.geom.index], lty=2)
  text(x=Tm.geom-5-273.15, y=17, labels="Tm", cex=2, col='black')
  segments(x0=Tm.geom-5-273.15, y0=16.5, x1=Tm.geom-273.15, y1=10)
  segments(x0=45, y0=0, x1=67, y1=0, lty=3, col='black')
  geom.line.fit <-  (dH.geom*(((((dH.geom/trim.data[[1]])-((dCp.geom*Tm.geom)/trim.data[[1]])+dCp.geom)
                                /(R*trim.data[[1]]))*K.geom)/((K.geom+1)^2))+((dCp.geom*K.geom)/(K.geom+1)))
  #lines(trim.data[[1]]-273.15, geom.line.fit,col= 'red')
  sigmoid <- c(rep(0, int.l), signal.sigmoid, rep(0, (length(trim.data[[1]]-T.end.index))))
  df <- data.frame(trim.data[[1]]-273.15, experiment, trim.data[[1]], sigmoid)
  colnames(df) <- c("Temp", "Cp.ave", "basetemp", "sigmoid")
  ggplot(df, aes(x=Temp, y = Cp.ave)) +
    geom_line(aes(x=Temp, y=Cp.ave)) +
    geom_line(aes(x=basetemp, y=sigmoid)) +
    geom_ribbon(data=subset(df, Temp < (T.end-273.15) & Temp > (t.upper-273.15) ) ,aes(ymin=sigmoid,ymax=Cp.ave), fill="blue", alpha="0.5")
  return(matrix(c("dH.raw", "dCp.raw", "Tm.raw", "dG(37).raw",
                  dH.geom,  dCp.geom,   Tm.geom,   dG37.geom), ncol= 4, byrow=TRUE))
}

############################################################################################
##########    End Graphical Features #######################################################
############################################################################################

############################################################################################
##########    Step 5. Van't Hoff raw data fit  #############################################
############################################################################################
Tm.geom <- 328.44
dH.VH              <- 4*R*(Tm.geom^2)*mean((diff(norm.int.min.sig[(Tm.geom-10):(Tm.geom+10)])/diff(trim.data[[1]][(Tm.geom-10):(Tm.geom+10)])))

###########################################################################################
########## Step 6. nonlinear regression modeling ##########################################
###########################################################################################
dH    <- as.numeric(Graphical.features(cp.ave)[2,1])           # give starting guess for the model using geometrical fit
dCp   <- as.numeric(Graphical.features(cp.ave)[2,2])
Tm    <- as.numeric(Graphical.features(cp.ave)[2,3])
HofT  <- dH+(dCp*(trim.data[[1]]-Tm))                          # define enthalpy using geom fit values
GofT  <- dH*(1-(trim.data[[1]]/Tm))+dCp*(trim.data[[1]]-Tm-trim.data[[1]]*log(trim.data[[1]]/Tm))
KofT  <- e^(-GofT/(R*trim.data[[1]]))
modfit <-  (dH*(((((dH/trim.data[[1]])-((dCp*Tm)/trim.data[[1]])+dCp)/(R*trim.data[[1]]))*KofT)/((KofT+1)^2))+((dCp*KofT)/(KofT+1)))
plot(trim.data[[1]], cp.ave, ylim=c(0,max(cp.ave)+1), pch=".") #PLOT data wo baseline to overlay graphical features to see how well it fits#########
lines(trim.data[[1]], modfit,col= 'red')                       # shows an estimate of the fit from the values guessed so far, Make sure its somewhat similar or the fitting will error out
par(mfrow=c(1,1), mar=c(4,7,2,1))

Temp <- rep(trim.data[[1]],3)
cp.ave <- c(data.wo.base[[1]],data.wo.base[[2]],data.wo.base[[3]])
df <- data.frame(Temp, cp.ave)
# Definition and nonlinear regression of the two state model fit. Model taken from Negative Coupling as a Mechanism for Signal Propagation between C2 Domains of Synaptotagmin I
# Michael E. Fealey,  Jacob W. Gauer,  Sarah C. Kempka,  Katie Miller,  Kamakshi Nayak,  R. Bryan Sutton,  Anne Hinderliter
two.state <- function(df) {
  dffit <- df
  fit.params <- nls(cp.ave~
                      (((dH+(dCp*(Temp-Tm)))*((dH/Temp)
                                              -((dCp*Tm)/Temp)+dCp)
                        *(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                              /(0.00198588*Temp))/(0.00198588*Temp))
                        /((exp(-((dH*(1-(Temp/Tm)))
                                 +(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))/(0.00198588*Temp))+1)^2))
                       +((dCp*(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                                   /(0.00198588*Temp))))
                         /(exp(-((dH*(1-(Temp/Tm)))+(dCp*(Temp-Tm-(Temp*log(Temp/Tm)))))
                               /(0.00198588*Temp))+1))), data=df,
                    start=list(Tm=Tm, dCp=dCp, dH=dH))
  sm.lm <- summary(fit.params)
  plot(dffit$Temp, dffit$cp.ave, xlab="Cp(T) (kcal/mol/K)", ylab="Temp(K)", pch=".")
  lines(dffit$Temp, predict(fit.params), col='orange')
  plot(dffit$Temp, dffit$cp.ave-predict(fit.params))
  dH    <- sm.lm$coefficients[3]
  dCp   <- sm.lm$coefficients[2]
  Tm    <- sm.lm$coefficients[1]
  dG37  <- (dH*(1-((273.15+37)/Tm))+dCp*((273.15+37)-Tm-(273.15+37)*log((273.15+37)/Tm)))
  regress.values <- list(sm.lm$parameters, dG37)
  dG <- dH*(1-(trim.data[[1]]/Tm))+dCp*(trim.data[[1]]-Tm-trim.data[[1]]*(log(trim.data[[1]]/Tm)))
  plot(trim.data[[1]], dG)
  confint2(fit.params, level=0.95, method=c("asymptotic"))
  conf <-nlsConfRegions(fit.params, length=300, exp=1)
  plot(conf, bounds=TRUE)
  cont <- nlsContourRSS(fit.params, lseq=20, exp=5)
  plot(cont, nlev=0, col=TRUE, col.pal=terrain.colors(100))
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
  par(mai = c(0.5,0.5, 0.5, 0.5))
  dconf <- data.frame(Tm = (conf$cr[,1]-273.15), dCp=conf$cr[,2], dH=conf$cr[,3], rss = conf$rss)
# Plotting confidence contours and intervals for pairwise and higher order evaluation of covariance of fit for goodness-of-fit
# evaluation.
  dconf %>%
    filter(rss < conf$rss95) %>%
    ggplot(aes(x=Tm, y=dCp))+
      geom_point(shape=4) +
      xlab(expression(paste("T"[m]," ( ",degree,"C)"))) +
      ylab(expression(paste(Delta,"C"[p]," (kcal/mol ",degree,"C)"))) +
      theme_bw()
  dconf%>% 
    filter(rss < conf$rss95) %>%
    ggplot(aes(x=Tm, y=dH))+
    geom_point(shape=4) +
    xlab(expression(paste("T"[m]," ( ",degree,"C)"))) +
    ylab(expression(paste(Delta,"H"[T[m]]," (kcal/mol)"))) +
    theme_bw()
  dconf %>%
    filter(rss < conf$rss95) %>%
  ggplot(aes(x=dCp, y=dH))+
    geom_point(shape=4) +
    xlab(expression(paste(Delta,"C"[p]," (kcal/mol ",degree,"C)"))) +
    ylab(expression(paste(Delta,"H"[T[m]]," (kcal/mol)"))) +
    theme_bw()
  return(fit.params)
}