
run_fish_coral <- function(time, env, pars) {
  
  # Get number of symbionts in run
  nsym <- length(pars$initS)  
  # Define synthesizing unit formula
  synth <- function(x, y, m) 1 / ((1 / m) + (1 / x) + (1 / y) - (1 / (x + y)))
  vsynth <- Vectorize(synth)  # Vectorize synthesizing unit for nsym > 1

  # Set initial values
  # ==================
  # Coral (host) fluxes
  ## Create empty vectors for each flux with length equal to time vector
  for(x in c("jX", "jN", "rNH", "rhoN", "jeC", "jCO2", "jHG", "jHT", "rCH", "dH.Hdt", "H")) {
    assign(x, rep(NA, length(time)))
  }
  ## Set initial values for host fluxes
  jX[1] <- (pars$jXm * env$X[1] / (env$X[1] + pars$KX))
  jN[1] <- (pars$jNm * env$N[1] / (env$N[1] + pars$KN))
  rNH[1] <- pars$jHT0 * pars$nNH * pars$sigmaNH
  rhoN[1] <- jN[1]
  jeC[1] <- 10
  jCO2[1] <- pars$kCO2 * jeC[1]
  jHG[1] <- 0.25
  jHT[1] <- pars$jHT0
  rCH[1] <- pars$jHT0 * pars$sigmaCH
  dH.Hdt[1] <- pars$jHGm
  H[1] <- pars$initH
  # Symbiont fluxes
  ## Create empty vectors for each flux with length equal to time vector
  for (x in c("rNS", "jL", "jCP", "jeL", "jNPQ", "jCO2w", "jSG", "rhoC", "jNw", "jST", "rCS", "cROS", "dS.Sdt", "S")) {
    assign(x, matrix(NA, ncol=nsym, nrow=length(time)))
  }
  ## Set initial values for symbiont fluxes
  rNS[1,] <- pars$jST0 * pars$nNS * pars$sigmaNS
  jL[1,] <- env$L[1] * pars$astar
  jCP[1,] <- pmax(0, vsynth(jL[1,] * pars$yCL, jCO2[1]*H[1]/pars$initS, pars$jCPm), na.rm=T)
  jeL[1,] <- pmax(jL[1,] - jCP[1,]/pars$yCL, 0)
  jNPQ[1,] <- pars$kNPQ
  jCO2w[1,] <- jCO2[1]*H[1]/pars$initS - jCP[1,]
  jSG[1,] <- pars$jSGm/10
  rhoC[1,] <- jCP[1,]
  jNw[1,] <- 0
  jST[1,] <- pars$jST0
  rCS[1,] <- pars$jST0 * pars$sigmaCS
  cROS[1,] <- 1
  dS.Sdt[1,] <- pars$jSGm
  S[1,] <- pars$initS
  # Initial fish and internal nitrogen conditions
  ## Create empty vectors for each flux with length equal to time vector
  for(x in c("dP.Pdt", "P", "dW.Wdt", "W", "dNi.dt", "Ni", "VH", "VHi", "M")) {
    assign(x, rep(NA, length(time)))
  }
  ## Set initial values of fish and internal nitrogen
  dP.Pdt[1] <- 0 
  P[1] <- pars$initP
  dW.Wdt[1] <- 0
  W[1] <- pars$initW
  VH[1] <- 1
  VHi[1] <- 0.5
  M[1] <- 0
  dNi.dt[1] <- 0
  Ni[1] <- env$N[1]

  # Run simulation by updating
  # ==========================
  dt <- time[2] - time[1]
  for (t in 2:length(time)) {

    # Photosynthesis fluxes
    # =====================
    # Light input flux
    jL[t,] <- (1.256307 + 1.385969 * exp(-6.479055 * (sum(S[t-1,])/H[t-1]))) * env$L[t] * pars$astar
    # CO2 input flux
    rCS[t,] <- pars$sigmaCS * (jST[t-1,] + (1-pars$yC)*jSG[t-1,]/pars$yC)  # metabolic CO2 recycled from symbiont biomass turnover
    # Production flux (photosynthetic carbon fixation)
    jCP[t,] <- vsynth(jL[t,] * pars$yCL, (jCO2[t-1] + rCH[t-1])*H[t-1]/sum(S[t-1,]) + rCS[t,], pars$jCPm) / cROS[t-1,]
    # Rejection flux: CO2 (wasted to the environment)
    jCO2w[t,] <- pmax((jCO2[t-1] + rCH[t-1])*H[t-1]/sum(S[t-1,]) + rCS[t,] - jCP[t,], 0)
    # Rejection flux: excess light energy not quenched by carbon fixation
    jeL[t,] <- pmax(jL[t,] - jCP[t,]/pars$yCL, 0)
    # Amount of excess light energy quenched by NPQ
    jNPQ[t,] <- (pars$kNPQ^(-1)+jeL[t,]^(-1))^(-1/1)  # single substrate SU
    # Scaled ROS production due to excess excitation energy (=not quenched by carbon fixation AND NPQ)
    cROS[t,] <- 1 + (pmax(jeL[t,] - jNPQ[t,], 0) / pars$kROS)^pars$k
    
    # Symbiont biomass fluxes
    # =======================
    # Nitrogen input flux
    rNS[t,] <- pars$jST0 * pars$nNS * pars$sigmaNS  # Recylced N from symbiont biomass turnover.
    # Production flux (symbiont biomass formation)
    jSG[t,] <- vsynth(pars$yC*jCP[t,], (rhoN[t-1]*H[t-1]/sum(S[t-1,]) + rNS[t,])/pars$nNS, pars$jSGm)
    # Rejection flux: carbon (surplus carbon shared with the host)
    rhoC[t,] <- pmax(jCP[t,] - jSG[t,]/pars$yC, 0)
    # Total amount of carbon shared by all symbionts
    rhoC.t <- sum(rhoC[t,]*S[t-1,])
    # Rejection flux: nitrogen (surplus nitrogen wasted to the environment)
    jNw[t,] <- max(rhoN[t-1]*H[t-1]/sum(S[t-1,]) + rNS[t,] - pars$nNS * jSG[t,], 0)
    # Total amount of nitrogen wasted by all symbionts
    jNw.t <- sum(jNw[t,]*S[t-1,])
    # Symbiont biomass loss (turnover)
    jST[t,] <- pars$jST0 * (1 + pars$b * (cROS[t,] - 1))
    
    # Host biomass fluxes
    # ===================
    # Food input flux (prey=both carbon and nitrogen)
    jX[t] <- (pars$jXm * env$X[t] / (env$X[t] + pars$KX))  # Prey uptake from the environment
    # Nitrogen input flux
    jN[t] <- (pars$jNm * Ni[t-1] / (Ni[t-1] + pars$KN))  # N uptake as a function of nitrogen concentration in coral head (Ni) 
    rNH[t] <- jHT[t-1] * pars$nNH * pars$sigmaNH  # Recycled N from host biomass turnover
    # Production flux (host biomass formation)
    jHG[t] <- synth(pars$yC*(rhoC.t/H[t-1] + jX[t]), (jN[t] + pars$nNX*jX[t] + rNH[t]) / pars$nNH, pars$jHGm)
    # Rejection flux: nitrogen (surplus nitrogen shared with the symbiont)
    rhoN[t] <- max(jN[t] + pars$nNX * jX[t] + rNH[t] - pars$nNH * jHG[t], 0)
    # Rejection flux: carbon -- given back to symbiont as CO2 input to photosynthesis
    jeC[t] <- max(jX[t] + rhoC.t/H[t-1] - jHG[t]/pars$yC, 0)
    # carbon not used in host biomass is used to activate CCM's that deliver CO2 to photosynthesis
    jCO2[t] <- pars$kCO2 * jeC[t] 
    # Host biomass loss
    jHT[t] <- pars$jHT0
    # metabolic CO2 recycled from host biomass turnover
    rCH[t] <- pars$sigmaCH * (jHT[t] + (1-pars$yC)*jHG[t]/pars$yC)

    # Convert coral (H) biomass to biovolume for interactions with fish / nitrogen
    VH[t] <- pars$kv * H[t-1]^pars$gamma
    VHi[t] <- VH[t] * pars$vi  # Internal volume of water
    # Coral mortality rate (breakage as a function of colony volume VH)
    M[t] <- pars$m * VH[t]^pars$mu

    # Balance equations
    # =================
    # Host (H)
    dH.Hdt[t] <- jHG[t] - pars$jHT0 - M[t] # Specific growth rates (Cmol/Cmol/d)
    # Symbiont (S)
    dS.Sdt[t,] <- jSG[t,] - jST[t,] - M[t] # Specific growth rates (Cmol/Cmol/d)
    # Damselfish (P)
    dP.Pdt[t] <- pars$rp * (pars$kp * VH[t] - pars$Bp * P[t-1] - pars$alpha.wp * W[t-1]) / (pars$kp * VH[t]) - pars$ap * env$U[t]
    # Hawkfish (W)
    dW.Wdt[t] <- pars$rw * (pars$kw * VH[t]^(2/3) - pars$Bw * W[t-1] - pars$alpha.pw * P[t-1]) / (pars$kw * VH[t]^(2/3)) - pars$aw * env$U[t]
    # Internal DIN concentration (Ni)
    dNi.dt[t] <- pars$D * (env$N[t] - Ni[t-1]) + (pars$ep*P[t-1] + pars$ew*W[t-1] + jNw.t*sum(S[t-1,]) - jN[t]*H[t-1])/VHi[t]

    # State variables
    # ===============
    H[t] <- H[t-1] + dH.Hdt[t] * H[t-1] * dt  # Biomass (Cmol)
    S[t,] <- S[t-1,] + dS.Sdt[t,] * S[t-1,] * dt  # Biomass (Cmol)
    P[t] <- P[t-1] + dP.Pdt[t] * P[t-1] * dt  # Biomass
    W[t] <- W[t-1] + dW.Wdt[t] * W[t-1] * dt  # Biomass
    Ni[t] <- Ni[t-1] + dNi.dt[t] * Ni[t-1] * dt # Ni concentration
  }

  # Return results
  # ==============
  out <- data.frame(
    time, env$L, env$N, env$X, env$U, jX=jX, jN=jN, rNH=rNH, rhoN=rhoN, jeC=jeC, 
    jCO2=jCO2, jHG=jHG, jHT=jHT, rCH=rCH, rNS=rNS, jL=jL, jCP=jCP, jeL=jeL, 
    jNPQ=jNPQ, jCO2w=jCO2w, jSG=jSG, rhoC=rhoC, jNw=jNw, jST=jST, rCS=rCS, 
    cROS=cROS, VH=VH, VHi=VHi, M=M, dH.Hdt=dH.Hdt, H=H, dS.Sdt=dS.Sdt, S=S, 
    dP.Pdt=dP.Pdt, P=P, dW.Wdt=dW.Wdt, W=W, dNi.dt=dNi.dt, Ni=Ni)
  return(out)
}
