rm( list=ls(all=TRUE) )


library(deSolve)
library(tidyr)
library(ggplot2)


###################################
##  setting
###################################


Tmin <- 0
Tmax <- 10.0

pars <- c(T0 = 10000, I0 = 100, V0 = 1000,
          lambda = 0, beta = 0.000002, d = 0,
          p = 10000, delta = 0.5, c = 10)



###################################
##  functions
###################################

system("R CMD SHLIB ode.c")
if (is.loaded("derivs")) {
  dyn.unload( paste("ode",.Platform$dynlib.ext,sep="") )
}
dyn.load( paste("ode",.Platform$dynlib.ext,sep="") )


ODEs <- function(pars, odemethod = 'rk4'){
  
  y <- c(Tp = as.numeric(pars["T0"]),
         Ip = as.numeric(pars["I0"]),
         Vp = as.numeric(pars["V0"]))
  
  times <- seq(Tmin, Tmax, step_size)
  out <- ode(y = y, parms = pars[-c(1:3)], times = times, func = "derivs", initfunc="initparms", nout = 1, outnames = c(""), dllname = "ode", odemethod = "rk4")
  
  as.data.frame(out)
}



###################################
##  fig
###################################

step_size <- 0.01
stime <- seq(Tmin,Tmax,step_size)

out <- list()

for (i in 1:2) {
  out[[i]] <- ODEs(pars = pars, odemethod = c("euler", "rk4")[i])
}

plt <- list()

for (i in 1:3) {
  
  plt[[i]] <- ggplot()
  
  for (j in 1:2) {
    
    data_temp <- cbind(out[[j]][, c("time", c("Tp", "Ip", "Vp")[i])], method = c("euler", "rk4")[j])
    colnames(data_temp) <- c("time", "value", "method")
    
    plt[[i]] <- plt[[i]] + 
      geom_line(data = data_temp, aes(x = time, y = value, linetype = method, color = method)) 
  }
  
  plt[[i]] <- plt[[i]] + 
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_color_manual(values = cbind(c("black", "red", "blue"), c("grey", "pink", "skyblue"))[i,]) +
    xlab("Time") + ylab("Value") + labs(title = c("Target cells", "Infected cells", "Viruses")[i]) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}


ggsave(wrap_plots(plt), filename = "plot_ode1.png", w = 12, h = 3)



step_size <- 1
stime <- seq(Tmin,Tmax,step_size)

out <- list()

for (i in 1:2) {
  out[[i]] <- ODEs(pars = pars, odemethod = c("euler", "rk4")[i])
}

plt <- list()

for (i in 1:3) {
  
  plt[[i]] <- ggplot()
  
  for (j in 1:2) {
    
    data_temp <- cbind(out[[j]][, c("time", c("Tp", "Ip", "Vp")[i])], method = c("euler", "rk4")[j])
    colnames(data_temp) <- c("time", "value", "method")
    
    plt[[i]] <- plt[[i]] + 
      geom_line(data = data_temp, aes(x = time, y = value, linetype = method, color = method)) 
  }
  
  plt[[i]] <- plt[[i]] + 
    scale_linetype_manual(values = c("solid", "dashed")) +
    scale_color_manual(values = cbind(c("black", "red", "blue"), c("grey", "pink", "skyblue"))[i,]) +
    xlab("Time") + ylab("Value") + labs(title = c("Target cells", "Infected cells", "Viruses")[i]) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}


ggsave(wrap_plots(plt), filename = "plot_ode2.png", w = 12, h = 3)

