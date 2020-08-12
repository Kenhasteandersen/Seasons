library(pracma)
source("plottools.R")
library(shiny)

ui <- fluidPage(
  
  # Application title
  titlePanel("Optimal seasonal strategies"),
  p("Implementation of the model in Thygesen et al: ",
    a(href="https://www.biorxiv.org/content/10.1101/2020.08.02.213595v1",
      'Optimal seasonal routines across animals and latitudes.')),
  p("Code on ",
    a(href="https://github.com/Kenhasteandersen/Seasons","github.")),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("rho",
                  "Seasonal amplitude:",
                  min = 0,
                  max = 2,
                  step = 0.05,
                  value = 1.5)
      ,
      sliderInput("w10",
                  "log10(body mass in g)",
                  min=-5,
                  max=5,
                  step = 0.25,
                  value=4)
    ),
    
    mainPanel(
      plotOutput("strategyPlot")
    )
  )
)

# Define server logic 
server <- function(input, output) {
  
  output$strategyPlot <- renderPlot({
    plotOptimal( calcOptimal(rho=input$rho, w=10^input$w10) )
  })
}

appSeasons = function() {
  shinyApp(ui = ui, server = server)
}


calcOptimal = function(rho, w) {
  # Time horizon for backward problem
  Tmax = 6  # Years
  dt = 0.001 # Years
  
  tvec = seq(0, Tmax, by=dt)
  nt = length(tvec)
  
  # State space (storage)
  Smax = 0.2
  nS = 31
  Svec = seq(0, Smax, length.out=nS)
  ds = Svec[2]-Svec[1]
  
  # Decision space (effort spent foraging)
  tauvec = seq(0,1,length.out = 51)
  ntau = length(tauvec)
  
  # Parameters
  
  # "Fixed parameters"
  fc0 = 0.1
  fc1 = 0.1
  a0 = 0.3
  a1 = 0.3
  R0 = 1.5
  
  # Initialize output
  #Vout = zeros(length(ws),length(rhos),length(Svec),length(tvec));
  #TauOut = Vout;
  #SigmaOut = Vout;
  #Stout = zeros(length(ws),length(rhos),length(tvec));
  #Vtout = Stout;
  #Tautout = Stout;
  #Sigmatout = Stout;
  
  
  wTCm = 1/40*w^0.25;
  
  # Available energy per effort as a function of time
  # Possible extension: Inlcuding an amplitude between 0 and 1, which depends on latitude
  R = function(t) R0*pmax(0,1-rho*cos(2*pi*t))
  
  # Available Energy as a function of foraging effort (functional response; also metabolism)
  E = function(t,tau) (tau*R(t))/(1+tau*R(t)) - fc0 - fc1*tau
  
  # Mortality as a function of foraging effort
  mu = function(t,tau) a0 + a1*tau
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Set up implementation
  
  # Grid fields
  V = matrix(data=0,nS,nt)
  I = V
  TAU = V
  SIGMA = V
  
  # Discretization of the gradient operator d/dS (right/up, finite differece)
  gradU = -diag(nS)/ds
  for (i in 1:nS-1) {
    gradU[i,i+1] = 1/ds
  }
  gradU[nS,nS] = 0
  
  # Discretization of the gradient operator d/dS (left/down, finite differece)
  gradD = diag(nS)/ds
  for (i in 1:nS-1) {
    gradD[i+1,i] = -1/ds
  }
  gradD[1,1] = 0
  
  # dV/dt + max{ dV/ds * (g-m) * sigma  - mu * V + (1-sigma)*g(-m) : tau, sigma}
  # sigma = 0: Reproduce
  # sigma = 1: Build storage
  # Negative available energy implies loosing storage
  
  # Backwards-in-time loop, solution of the DP equation
  for (i in seq(nt,2,by=-1)) {
    # Dirichlet boundary condition
    V[1,i] = 0
    
    # Rate of change in fitness assuming all surplus energy is allocated to
    # reproduction
    dV0 = matrix(rep(1,nS)) %*% t(matrix( pmax(E(tvec[i-1] ,tauvec),0) ))            # Reproducing
    dV0 = dV0 + (gradD%*%V[,i]) %*% t(matrix(pmin(E(tvec[i-1],tauvec) ,0)))    # Withdrawing from storage
    dV0 = dV0 - V[,i]  %*% t(matrix(mu(tvec[i-1],tauvec)))                     # Mortality
    
    # Rate of change in fitness assuming all surplus energy is allocated to
    # building storage
    dV1 = (gradU%*%V[,i]) %*% t(matrix(pmax(E(tvec[i-1],tauvec),0)))            # Building storage
    dV1 = dV1 + (gradD%*%V[,i]) %*% t(matrix( pmin(E(tvec[i-1],tauvec) ,0)))   # Withdrawing from storage
    dV1 = dV1 - V[,i] %*% t(matrix(mu(tvec[i-1],tauvec)))                    # Mortality
    
    # Choose optimal storage allocation, 0 or 1
    # TODO: Allow the critter to follow the switching surface in stead of jittering back and forth over it
    dV = pmax(dV0,dV1);
    
    # Choose optimal effort level
    I[,i-1] = max.col(dV, ties.method = "first")
    dV = apply(dV,1,max)
    V[,i-1] = V[,i] + dV*dt/wTCm;
    V[1,i-1] = 0
    
    TAU[,i-1] = tauvec[I[,i-1]];
    
    # Now find again the optimal allocation, for the optimal effort level
    # TODO: Vectorize this!
    for (j in 1:nS)
      SIGMA[j,i-1] = (max(dV1[j,])>=max(dV0[j,]))
    
    # disp(max(V(:,i-1)))
  }
  
  # Save output
  Vout = V;
  TauOut = TAU;
  SigmaOut = SIGMA;
  
  
  
  # if DoPlot
  #     figure(1)
  #     imagesc(tvec,Svec,V)
  #     set(gca,'ydir','normal')
  #     title('Fitness')
  #     colorbar
  #     
  #     figure(2)
  #     imagesc(tvec,Svec,TAU)
  #     set(gca,'ydir','normal')
  #     title('Effort')
  #     colorbar
  #     
  #     figure(3)
  #     imagesc(tvec,Svec,SIGMA)
  #     set(gca,'ydir','normal')
  #     title('Storage  building')
  #     colorbar
  #     
  #     drawnow
  # end
  
  
  # Solve forward in time
  
  # Initialize solutions
  St = vector(length=length(tvec), mode="numeric")   # State (storage)
  Et = St;                  # Available energy
  Sigmat = St;              # Chosen allocation
  Taut = St;                # Chosen effort
  mut = St;                 # Mortality
  Vt = St;                  # Fitness
  
  St[1] = Svec[2]
  
  for (i in 1:length(tvec)) {
    Taut[i] = interp1(Svec,TAU[,i],St[i],method='nearest');
    Sigmat[i] = interp1(Svec,SIGMA[,i],St[i],method='nearest');
    Et[i] = E(tvec[i],Taut[i]);
    mut[i] = mu(tvec[i],Taut[i]);
    Vt[i] = interp1(Svec,V[,i],St[i],method='nearest');
    
    if (i<length(tvec))
      St[i+1] = max(0,St[i] + Et[i]*Sigmat[i] * dt / wTCm);
    
  }
  
  dSt = diff(St)/dt;
  Resource = R(tvec)
  return( list(t=tvec, tau=Taut, sigma=Sigmat, V=Vt, S=St, E=Et, R=Resource) )
}

plotOptimal = function(opt) {
  ix = opt$t>4 & opt$t<5
  t = opt$t[ix]-4
  
  defaultplotvertical( nPanels=2 )
  
  defaultpanel(xlim=c(0,1), ylim=c(0,4), ylab="Resource", xaxis = FALSE)
  ribbon(t, ymax=opt$R[ix], ymin=0*t, col="blue")
  lines(t,rep(1.5, length(t)), lty=dashed)
  
  
  defaultpanel(xlim=c(0,1), ylim=c(0,1),
               xlab="Time", ylab="Fraction")
  
  hibernation = (opt$tau==0) & (opt$V>0)
  ribbon(t, ymax=hibernation[ix], ymin=0*t,col=stdgrey)

  act = opt$tau[ix]
  lines(t, act, col="red", type="l", lwd=2)
  
  sigma = smoothing(opt$sigma[ix])
  sigma[act==0] = NA
  lines(t, sigma, col="blue", lwd=2)
  ribbon(t, ymax=opt$S[ix], ymin=0*t,col="blue")
  
  dying = opt$V[ix]<=0
  ribbon(t, ymax=dying, ymin=0*t,col=darkgrey)
  #lines(t, opt$S[ix], col="magenta", lwd=3)
  #lines(t, opt$E[ix], col="brown", lwd=3)
  
  legend("bottom", 
         legend=c("Fitness<0","Hibernation","Storage",
                  "Activity","Accumulate storage"),
         col=c(NA, NA, NA,"red", "blue"),
         lty=c(0,0,0,1,1),
         lwd=2,
         bty="n",
         fill=c(darkgrey, stdgrey, "blue", "transparent", "transparent"),
         border=c("transparent")
         
  )
}

smoothing = function(y) {
  lgt = 60; # Width of window used to search of flutuating solution in sigma
  windowSize = 40; # Width of window used for smoothing
  
  sigmaForward = 0*y;
  sigmaBackward = 0*y;
  for (i in 1:length(y)) {
    sigmaBackward[i] = mean(y[max(1,i-lgt):i]);
    sigmaForward[i] = mean(y[i:min(length(y),i+lgt)]);
  }
  ixSmooth = ((sigmaForward != 0) & (sigmaBackward != 0)) + 
    ((sigmaForward != 1) & (sigmaBackward !=1))
  
  ii = 1
  while( first(ixSmooth[ii:length(ixSmooth)],2)>=ii) {
    iStart = first(ixSmooth[ii:length(ixSmooth)],2)
    iEnd = first(ixSmooth[iStart:length(ixSmooth)],1)+iStart
    
    window = min(windowSize, iEnd-iStart)
    y[iStart:iEnd] = filter(y[iStart:iEnd], rep(1 / window, window), 
                            sides = 2, circular=TRUE)
    
    #   filter((1/windowSize)*ones(1,windowSize), 1, y[iStart:iEnd]);
    ii = iEnd;
  }
  return(y)
  
}

first = function(x,y) {
  z = which(x==y)
  z = z[1]
  if (is.na(z))
    z = FALSE
  return(z)
}
