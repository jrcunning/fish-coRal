
def_pars <- function(nsym=1,
                     jHT0=0.03,  # Host specific biomass turnover rate (d^-1)
                     nNH=0.18,  # N:C ratio in host biomass (-)
                     nNX=0.2,  # N:C ratio in prey biomass (-)
                     sigmaNH=0.9,  # Proportion of host nitrogen turnover recycled (-)
                     sigmaCH=0.1,  # Proportion of host carbon turnover recycled (-)
                     jXm=0.13,  # Maximum specific host feeding rate (molX/CmolH/d)
                     jNm=0.035,  # Maximum specific host DIN uptake rate (molN/CmolH/d)
                     jHGm=1,  # Maximum specific host growth rate (CmolH/CmolH/d)
                     kCO2=10,  # Rate of host CCM's (molCO2/molC/d)
                     KN=1.5e-6,  # Half-saturation constant for host DIN uptake (molN/L)
                     KX=1e-6,  # Half-saturation constant for host feeding (CmolX/L)
                     initH=1,  # Initial host biomass (CmolH)
                     yC=0.8,
                     jST0=rep(0.03, nsym),  # Symbiont specific biomass turnover rate (d^-1)
                     nNS=rep(0.13, nsym),  # N:C ratio in symbiont biomass (-)
                     yCL=rep(0.1, nsym),  # L:C ratio in fixed carbon (=quantum yield) (molC/mol ph)
                     kNPQ=rep(112, nsym),  # capacity of non-photochemical quenching (mol ph/CmolS/d)
                     # calculated as 4x max. photochemical quenching (Gorbunov et al. 2001)
                     kROS=rep(80, nsym),  # amount of excess light beyond NPQ capacity (e.g., jeL-jNPQ) that doubles ROS production relative to baseline (mol ph/CmolS/d)
                     k=rep(1, nsym),  # exponent on ROS production (-)
                     astar=rep(1.34, nsym),  # Symbiont specific cross-sectional area (m^2/C-molS)
                     sigmaNS=rep(0.9, nsym),  # Proportion of symbiont nitrogen turnover recylced (-)
                     sigmaCS=rep(0.9, nsym),  # Proportion of symbiont carbon turnover recycled (-)
                     jCPm=rep(2.8, nsym),  # Maximum specific photosynthate production rate (Cmol/CmolS/d)
                     jSGm=rep(0.25, nsym),  # Maximum specific symbiont growth rate (CmolS/CmolS/d)
                     initS=rep(initH/10, nsym),  # Initial symbiont biomass (CmolS)
                     b=rep(5, nsym),  # Scaling parameter for bleaching response
                     #
                     kv=3.068,  # liters per C-mol H
                     gamma=1.603,
                     m=0.002, # 1/(d * (cubic meters)^mu)
                     mu=1,
                     # Damselfish parameters (P)
                     kp=1,  # carrying capacity
                     initP= kv * initH ^ gamma * kp, # (kv * initH ^ gamma) * kp,  # initial biomass
                     rp=.05, # intrinsic growth rate
                     Bp=1,
                     alpha.wp=1,
                     ap=1,  # attack rate on damselfish
                     # Hawkfish parameters (W)
                     kw=1,
                     initW=(kv * initH ^ gamma)^2/3 * kw,  # (kv * initH ^ gamma)^2/3 * kw,  # initial biomass
                     rw=.05, # intrinsic growth rate
                     Bw=1,
                     alpha.pw=1,
                     aw=1,  # attack rate on hawkfish
                     U=1,  # Uber-predator abundance
                     ep=.001,  # damselfish N excretion rate (N per biomass)
                     ew=.001,  # hawkfish N excretion rate (N per biomass)
                     D=360,  # flushing rate of N between interstitial volume and ambient environment (per day)
                     vi=0.5) {  # fraction of coral volume occupied by water) {
  return(list(
    jHT0=jHT0, nNH=nNH, nNX=nNX, sigmaNH=sigmaNH, sigmaCH=sigmaCH, jXm=jXm, jNm=jNm, jHGm=jHGm, 
    kCO2=kCO2, KN=KN, KX=KX, initH=initH, yC=yC, jST0=jST0, nNS=nNS,
    yCL=yCL, kNPQ=kNPQ, kROS=kROS, k=k, astar=astar, sigmaNS=sigmaNS, sigmaCS=sigmaCS, jCPm=jCPm, 
    jSGm=jSGm, initS=initS, b=b, kv=kv, gamma=gamma, m=m, mu=mu, kp=kp,
    initP=initP, rp=rp, Bp=Bp, alpha.wp=alpha.wp, ap=ap, kw=kw, initW=initW, rw=rw, Bw=Bw, 
    alpha.pw=alpha.pw, aw=aw, U=U, ep=ep, ew=ew, D=D, vi=vi
  ))
}

