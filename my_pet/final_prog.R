Full Bioaccumulation Program:
  
  ############
# PACKAGES #
############
# Package solving differential equations in R
# install.packages("deSolve")   
library(deSolve)



########################
# DATA EXP 1 SOLEBEMOL #
########################

# 6 months old fish + 84 days contamination + 84 days elimination
# exp="1PCB"
# T_C=19 and T_N=NA => constant temperature
# n_d1=6*30  => 6 months old fish at t0 of the experiment
# n_d2=84   => 84 days contamination
# n_d3=84   => 84 days elimination
# rho_food_wet=21612 => pellets

# Food measurements for PCB modality:
exp_deconta3_PCB_food <- read.delim2("exp_Solebemol1(3+3deconta)_PCB_food_RFlo.txt",h=T, sep="\t", dec=",", na.strings = "NA")
# as the table contains 1 value per day we have to create a value for each time step of integration with a time step of 0.5 days
# duplication of each value
deconta3_PCB_food <- seq(from=0, to=exp_deconta3_PCB_food[1,2], length.out=6*30*2+1) # linear augmentation of food from 0 to the firste measured value
for(i in 1:dim(exp_deconta3_PCB_food)[1]){
  deconta3_PCB_food <- c(deconta3_PCB_food, replicate(exp_deconta3_PCB_food[i,2],n=2) )
}
plot(deconta3_PCB_food)

# Biometry for all modalities:
exp_deconta3_biometry <- read.delim2("exp_Solebemol1(3+3deconta)_biometry_all_modilities_RFlo.txt",h=T, sep="\t", dec=",", na.strings = "NA")
exp_deconta3_PCB_biometry <- exp_deconta3_biometry[exp_deconta3_biometry$modality=="PCB",]
exp_deconta3_CS_biometry <- exp_deconta3_biometry[exp_deconta3_biometry$modality=="control_solvant",]

# Levels of contamination PCB:
exp_deconta3_PCB_contam <- read.delim2("exp_Solebemol1(3+3deconta)_PCB_contam_RFlo.txt",h=T, sep="\t", dec=",", na.strings = "NA")
exp_deconta3_PCB_contam$CB105<-as.numeric(as.character(exp_deconta3_PCB_contam$CB105))
exp_deconta3_PCB_contam$CB118<-as.numeric(as.character(exp_deconta3_PCB_contam$CB118))
exp_deconta3_PCB_contam$CB149<-as.numeric(as.character(exp_deconta3_PCB_contam$CB149))
exp_deconta3_PCB_contam$CB153<-as.numeric(as.character(exp_deconta3_PCB_contam$CB153))


########################
# DATA EXP 2 SOLEBEMOL #
########################


## TEMPERATURE ##
exp2_TS_temperature <- read.delim2("exp2_TS_temperature_RFlo.txt", h=T, sep="\t", dec=",")
# as the table contains 1 value per MONTH we have to create a value for each time step of integration with a time step of 0.5 days
# *30*2 of each value
exp2_temperature <- replicate(exp2_TS_temperature$temperature[1],n=30/0.5)
for(i in 2:dim(exp2_TS_temperature)[1]){
  exp2_temperature <- c(exp2_temperature, replicate(exp2_TS_temperature$temperature[i],n=30/0.5) )
}
plot(exp2_temperature, type="l",
     main="Time evolution of temperature during the experiment",
     xlab="Time since fecundation (days)", ylab="Experimental temperature (°C)")


## FOOD ##
exp2_TS_food <- read.delim2("exp2_TS_food_RFlo.txt", h=T, sep="\t", dec=",")
# as the table contains 1 value per MONTH we have to create a value for each time step of integration with a time step of 0.5 days
# *30*2 of each value
exp2_food <- seq(from=0, to=exp2_TS_food$joules_ingérés_modif[1], length.out=6*30*2+1) # linear augmentation of food from 0 to the firste measured value
for(i in 1:dim(exp2_TS_food)[1]){
  exp2_food <- c(exp2_food, replicate(exp2_TS_food$joules_ingérés_modif[i],n=30/0.5) )
}
plot(exp2_food, type="l",
     main="Time evolution of ingestion during the experiment (data modif)",
     xlab="Time since fecundation (days)", ylab="Ingestion (Joules per day per fish)")


## BIOMETRY ##



## PCB CONTAMINATION ##
exp2_PCB_contamination <- read.delim2("exp2_PCB_contamination_RFlo.txt", h=T, sep="\t", dec=",")


# ELIM
# 6 months old fish + 12 months contamination + 6 months elimination
# exp="PCBelim"
# T_C=0 and T_N=exp2_temperature => measured temperature
# n_d1=6*30    => 6 months old fish at t0 of the experiment
# n_d2=12*30   => 12 months contamination
# n_d3=6*30    => 6 months elimination
# rho_food_wet=21612 => pellets
exp2_PCB_elim <- subset(exp2_PCB_contamination,Phase=="conta")
exp2_PCB_elim <- rbind(exp2_PCB_contamination[exp2_PCB_contamination$Phase=="elim",],
                       exp2_PCB_elim[exp2_PCB_elim$time_d<450,])

# CONTA18
# 6 months old fish + 18 months contamination
# exp="PCBconta18"
# T_C=0 and T_N=exp2_temperature => measured temperature
# n_d1=6*30    => 6 months old fish at t0 of the experiment
# n_d2=18*30   => 18 months contamination
# n_d3=0    => no elimination
# rho_food_wet=21612 => pellets
exp2_PCB_conta <- subset(exp2_PCB_contamination,Phase=="conta")
exp2_PCB_conta <- subset(exp2_PCB_conta,time_d<570)


# CONTA6y
# 6 months old fish + 6 years contamination
# exp="PCBconta6y"
# T_C=0 and T_N=exp2_temperature => measured temperature
# n_d1=6*30    => 6 months old fish at t0 of the experiment
# n_d2=6*12*30   => 6 years contamination
# n_d3=0    => no elimination
# rho_food_wet=21612 => pellets
exp2_PCB_conta6y <- subset(exp2_PCB_contamination,Phase=="conta")




##################
# BIOACC PROGRAM #
##################



##########################
## FIXED DEB PARAMETERS ##
##########################

## AUXILIARY PARAMETERS ##

deltaM=0.204         # shape coefficient (-)
dV=1                 # density of structure (gw / cm^3)

C_V_w=0.814          # water content of structure (-)
C_E_w=0.725          # water content of reserve (-)
C_Gm_w=0.85          # water content of male gonads (-)
C_Gf_w=0.66          # water content of female gonads (-)

rho_V_dry=20070      # energy to dry weight of structure conversion (J/gd)
rho_E_dry=25940      # energy to dry weight of reserve conversion (J/gd)
rho_Gm_dry=20000     # energy to dry weight of male gonads conversion (J/gd)
rho_Gf_dry=20720     # energy to dry weight of female gonads conversion (J/gd)

Negg_per_gw=519      # number of eggs in 1g of female total wet weigth (#/gw)

W_E_dry_0=0.00012    # dry weight of reserve in the egg (gd)
Lw_0=0.14            # initial physical length (cm)
Lw_pm=22             # physical length at maturity for male (cm)
Lw_pf=31             # physical length at maturity for female (cm)


## PRIMARY PARAMETERS ##

Tref=273+10          # reference temperture for fluxes parameters (K)
TA=4550              # Arrhenius temperature (K)
pXm=460              # maximum ingestion rate at Tref (J.cm^-2.d^-1)
kappa_X=0.8          # assimilation efficiency (-)
pM=18.1              # volume-specific somatic maintenance cost at Tref (J.cm^-3.d^-1) 
EG=7000              # volume-specific cost of structure (J.cm^-3)
Em=2903              # maximum reserve density (J.cm^-3)
kappa_m=0.64         # fraction of mobilised reserve allocated to soma for male (-)
kappa_f=0.7          # fraction of mobilised reserve allocated to soma for female (-)
kappa_R_m=0.001      # fraction of reproduction energy fixed in male gamets (-)
kappa_R_f=0.5        # fraction of reproduction energy fixed in female gamets (-)

## COMPOUND PARAMETERS ##

Vp_m = (Lw_pm*deltaM)^3    # volume at first maturity for male (cm^3)
Vp_f = (Lw_pf*deltaM)^3    # volume at first maturity for female (cm^3)

## STATE VARIABLES AND INITIAL VALUES ##

E_0 = W_E_dry_0 * rho_E_dry   # energy in the reserve compartment at the beginning of the development (J)

V_0 = (Lw_0 * deltaM)^3       # structural volume at the beginning of the development (cm3)

GAM_cum_0 = 0                  # energy in the reproduction compartment at the beginning of the development (J)

Q_Wcum_0 = 0                  # PCB congener content at the beginning of the development (ng)

inits = c(E_0,V_0,GAM_cum_0, Q_Wcum_0)   # initial conditions vector for variables integrated


#################################
## TIME VECTOR FOR INTEGRATION ##
#################################

a_0=0.5                           # Age since fecondation at the beginning of the simulation
t_step=0.5                        # Time step





BIOACC <- function(exp,              # identification of the experiment
                   contaminant,      # name of the contaminant
                   n_d1,                # number of days before the experiment
         n_d2,			 # number of days of contamination (d)
         n_d3,			 # number of days of decontamination (d)
         rho_food_wet, 	# energy to wet weight of food conversion
         T_C,		# constant temperature (°C) => 0 if non constant T°
         T_N,            # vector containing temperature kinetic (°C) => NA if constant T°
         ae, 		# assimilation efficiency of the contaminant (-)
         log_Kow,        # octanol-water partition coefficient (-)
         C_X,		# concentration of PCB in food (ng.gd)
         C_X_res, 	# residual contamination of non contaminated food (ng.gd)
         J_food,    # vector containing food kinetic (J/d.fish) => NA if constant food
         food_f,    # 
         food_m,
         ingestion)		# type of information used for ingestion: 0 for f, 1 for joules per day per fish
{
  
  te = 3 + 0.033 * log_Kow       # transfer efficiency of PCB from body to eggs (-)
  t=seq(from=a_0, to= a_0 + n_d1 + n_d2 + n_d3, by=t_step)   # integration time vector
  
  parameters <- list(deltaM=deltaM, dV=dV, C_V_w=C_V_w, C_E_w=C_E_w, C_Gm_w=C_Gm_w, C_Gf_w=C_Gf_w,
                     rho_food_wet=rho_food_wet,
                     n_d1=n_d1,n_d2=n_d2,n_d3=n_d3,
                     rho_V_dry=rho_V_dry, rho_E_dry=rho_E_dry, rho_Gm_dry=rho_Gm_dry,
                     rho_Gf_dry=rho_Gf_dry, Negg_per_gw=Negg_per_gw, W_E_dry_0=W_E_dry_0,
                     Lw_0=Lw_0, Lw_pm=Lw_pm, Lw_pf=Lw_pf, Tref=Tref,TA=TA,
                     pXm=pXm, kappa_X=kappa_X, pM=pM, EG=EG, Em=Em,kappa_m=kappa_m, 
                     kappa_f=kappa_f,kappa_R_m=kappa_R_m,kappa_R_f=kappa_R_f,
                     J_food=J_food, food_f=food_f, food_m=food_m, Vp_m=Vp_m, Vp_f=Vp_f,
                     ae=ae, te=te, T_C=T_C, T_N=T_N,
                     E_0=E_0, V_0=V_0, GAM_cum_0=GAM_cum_0,
                     Q_Wcum_0=Q_Wcum_0,
                     a_0=a_0, t_step=t_step,
                     C_X=C_X, C_X_res=C_X_res, log_Kow=log_Kow)
  
  ##########################
  ## INTEGRATION FUNCTION ##
  ##########################
  
  integration <- function(t,f,parms){
    
    
    ## FORCING VARIABLES ##                          
    
    CX = c(rep(C_X_res,times=1+n_d1/t_step), rep(C_X,times=n_d2/t_step), rep(C_X_res,times=n_d2/t_step) )
    
    if(T_C==0){
      T = 273 + T_N[t/t_step]
    } else {T = 273+T_C}
    
    
    ## TEMPERATURE CORRECTION ##
    
    c_T = exp(TA/Tref - TA/T)
    
    pXm_T = pXm * c_T     # T° corrected maximum ingestion rate
    IXm_T = pXm_T / rho_food_wet   # T° corrected maximum ingestion rate in grams
    pM_T = pM * c_T       # T° corrected volume-specific somatic maintenance cost
    
    
    ## FLUXES ##
    if(ingestion==1){
      pX = J_food[t/t_step]                               # ingestion food disp non cte
    } else pX = pXm_T * parms["food"] * f[2]^(2/3)                    # ingestion food disp cte  

    pA = kappa_X * pX                                            # assimilation
    pC = ( (EG/Em) * kappa_X * pXm_T * f[2] ^(-1/3) + pM_T) / 
      ( EG/f[1] + parms["kappa"]/f[2] )                          # catabolization
    pM = pM_T * f[2]                                             # structural maintenance
    pG = parms["kappa"] * pC - pM                                # growth
    pJ = ((1 - parms["kappa"]) / parms["kappa"]) * pM_T * min (f[2],parms["Vp"])   # maturity maintenance
    pR = (1-parms["kappa"]) * pC - pJ                            # maturity / reproduction
    
    
    ## DIFFERENTIAL EQUATIONS ##
    
    df1 = pA - pC                                       # f1 = reserve (E)
    df2 = pG / EG                                       # f2 = structure (V)
    df3 = parms["kappa_R"] * pR * (f[2]>=parms["Vp"])   # f3 = reproduction (GAM)
    
    if(ingestion==1){
      df4 = pX/rho_food_wet * ae * CX[t/t_step]
    }else df4 = IXm_T * parms["food"] * f[2]^(2/3) * ae * CX[t/t_step] # f4 = assimilated PCB (Q_Wcum)
    
    
    
    ## OUTPUTS OF THE INTEGRATION FUNCTION ##
    
    list(c(df1,df2,df3,df4))                 # output of the function
    
  }
  
  ###############
  ## SOLUTIONS ##
  ###############
  
  states <- function(food, kappa, Vp, kappa_R){
    
    solutions <- lsoda(y=inits, times=t , func=integration, 
                       parms=c(food=food,kappa=kappa,Vp=Vp, kappa_R=kappa_R))
    E <- solutions[,2]         # time evolution of reserve (E, in J)
    V <- solutions[,3]         # time evolution of structural volume (V, in cm3)
    GAM_cum <- solutions[,4]    # time evolution of cumulated energy for reproduction (GAM, in J)
    Q_Wcum <- solutions[,5]    # time evolution of cumulated assimilated PCB (Q_Wcum, in ng of PCB)
    
    # time evolution of energy contained within the gamets (J)
    GAM <- rep.int(NA,times=length(t))         # creation of empty vector
    GAM[1] <- GAM_cum_0                         # initial value = 0
    
    for(i in 2:length(t)){    # starting from i=2 permits to use [i-1] as a position in the vector
      # The day of spawing :
      if(t[i]%%365==0){GAM[i] = 0}    # spawning = state variable of energy in gamets is artificialy emptied
      # The rest of the year :
      else {GAM[i] = GAM[i-1] + GAM_cum[i] - GAM_cum[i-1]}} # state variable of energy in gamets
    # = energy in gamets just before
    # + energy allocated to gamets between the 2 time steps
    
    # Energy contained within gamets just before reproduction (J)
    GAM_prior_repro <- rep.int(NA,times=length(t))
    for(i in 1:length(t)){
      if(t[i]%%364==0) {GAM_prior_repro[i] <- GAM[i]}
      else GAM_prior_repro[i] <- NA}
    
    list(E=E, V=V, GAM_cum=GAM_cum, Q_Wcum=Q_Wcum, GAM=GAM, GAM_prior_repro=GAM_prior_repro)
    
  }
  
  state_m <- states(food=food_m, kappa=kappa_m,Vp=Vp_m, kappa_R=kappa_R_m)
  state_f <- states(food=food_f, kappa=kappa_f,Vp=Vp_f, kappa_R=kappa_R_f)
  
  
  ##########################
  ## OBSERVABLE VARIABLES ##
  ##########################
  
  obs <- function(E, V, Q_Wcum, GAM, GAM_prior_repro, kappa, Vp, rho_G_dry, C_G_w){
    
    # time evolution of physical length (Lw, in cm)
    Lw <- V^(1/3) / deltaM
    
    # time evolution of wet weight of structure (in gw)
    WVw <- dV * V
    
    # time evolution of dry weight of structure (in gd)
    WVd <- WVw * (1-C_V_w)
    
    # time evolution of wet weight of reserve (in gw)
    WEw <- E / (rho_E_dry * (1-C_E_w))
    
    # time evolution of dry weight of reserve (in gd)
    WEd <- E / rho_E_dry
    
    # dry weight of gamets (gd)
    WGAMd = GAM / rho_G_dry
    
    # wet weight of gamets (gw)
    WGAMw = GAM / (rho_G_dry * (1-C_G_w))
    
    # Gamets dry weight just before reproduction (gd)
    WGAMd_prior_repro <- GAM_prior_repro / rho_G_dry
    
    # Gamets wet weight just before reproduction (gw)
    WGAMw_prior_repro <- GAM_prior_repro / (rho_G_dry * (1-C_G_w))
    
    # total dry weight (gd)
    WWd = WEd + WVd + WGAMd
    
    # total wet weight as sum of wet weight of compartments (gw)
    WWw_sum = WEw + WVw + WGAMw
    
    # Fulton's dry index:
    K_dry <- 100*WWd / (Lw)^3
    
    # Estimated dry content from Fonds et al. 1989 (-)
    dry_content <- 40.68 * K_dry^0.364 / 100
    
    # total wet weight from dry weight corrected by dry content from Fulton's relation (gw)
    WWw_Fulton <- WWd / dry_content
    
    # Fulton's wet index:
    Kw_WWw_sum <- 100*WWw_sum/Lw^3
    Kw_WWw_Fulton <- 100*WWw_Fulton/Lw^3
    
    # GAMado-somatic index (GSI) :
    GSI_WWw_sum = 100 * WGAMw / WWw_sum
    GSI_WWw_Fulton = 100 * WGAMw / WWw_Fulton
    
    ########################
    ## Energy in total body
    
    # Calculation of EW_per_gd from simulations:
    EW_per_gd <- rep.int(NA,times=length(t))
    for(i in 1:length(t)){ #from 2 to avoid t=0 and division by 0
      if(V[i]<Vp) EW_per_gd[i] = (E[i] + WVd[i]*rho_V_dry) / WWd[i]
      else EW_per_gd[i]=NA}
    
    # Calculation of EW_per_gw from simulations :
    EW_per_gw_WWw_sum <- rep.int(NA,times=length(t))
    for(i in 1:length(t)){ #from 2 to avoid t=0 and division by 0
      if(V[i]<Vp) EW_per_gw_WWw_sum[i] = (E[i] + WVd[i]*rho_V_dry) / WWw_sum[i]
      else EW_per_gw_WWw_sum[i]=NA}
    
    EW_per_gw_WWw_Fulton <- rep.int(NA,times=length(t))
    for(i in 1:length(t)){ #from 2 to avoid t=0 and division by 0
      if(V[i]<Vp) EW_per_gw_WWw_Fulton[i] = (E[i] + WVd[i]*rho_V_dry) / WWw_Fulton[i]
      else EW_per_gw_WWw_Fulton[i]=NA}
    
    
    ###############################
    ## Energy density of structure
    
    # Calculation of rho_V_dry with a given value of EW_per_gw for juvenile:
    # If we define the total body energy density as a function of Fulton's dry index (EW_per_gd_value = 10^3*35.69*Kd_m^0.273, Fonds 1989) :
    rho_V_dry_value <- rep.int(NA,times=length(t))
    for(i in 1:length(t)){
      if(V[i]<Vp) rho_V_dry_value[i] = (10^3*35.69*K_dry[i]^0.273 * WWd[i] - E[i]) / WVd[i]
      else rho_V_dry_value[i]=NA}
    
    
    ###########################################
    ## Calculation of total wet concentration:
    
    humidity_sum <- 1-(WWd / WWw_sum)
    humidity_Fulton <- 1-(WWd / WWw_Fulton)
    
    
    ###########################################
    ## Reserve proportion :
    EonWd <- WEd/WWd
    EonWw_sum <- WEw/WWw_sum
    EonWw_Fulton <- WEw/WWw_Fulton
    
    #########################
    ## SPAWNING FRO FEMALE ##
    
    # From energy of the gamets buffer:
    
    Negg_buffer_energy <- GAM_prior_repro / E_0
    
    Negg_buffer <- WGAMd_prior_repro / W_E_dry_0
    
    # From the fecundity value of Deniel with total wet weight:
    # Creation of empty vector :
    Negg_sum <- rep.int(NA,times=length(t))         # creation of the vector of number of eggs layed for the actual year
    Negg_sum[1] <- 0                                # initial value = 0
    Negg_Fulton <- rep.int(NA,times=length(t))         # creation of the vector of number of eggs layed for the actual year
    Negg_Fulton[1] <- 0                                # initial value = 0
    
    for(i in 2:length(t)){    # starting from i=2 permits to use [i-1] as a position in the vector
      
      # The day of spawing :
      if(t[i]%%365==0){
        if(V[i]<Vp_f) {Negg_sum[i] = 0
                       Negg_Fulton[i] = 0}
        else {Negg_sum[i] = WWw_sum[i-1]*Negg_per_gw     # number of eggs layed for the actual year
              Negg_Fulton[i] = WWw_Fulton[i-1]*Negg_per_gw}
      }
      # The rest of the year :
      else {Negg_sum[i] = NA                           # no spawning
            Negg_Fulton[i] = NA
      }
    }
    
    #########################
    ## BIOACCUMULATION ##
    
    Q_GAM <- rep.int(NA,times=length(t))       # creation of the vector of 
    Q_GAM[1] <- 0                              # initial value
    
    C_GAM <- rep.int(NA,times=length(t))       # creation of the vector of 
    C_GAM[1] <- 0                              # initial value
    
    Q_W <- rep.int(NA,times=length(t))         # creation of the vector of 
    Q_W[1] <- Q_Wcum_0                         # initial value
    
    for(i in 2:length(t)){ 
      Q_GAM[i] = 0 + (V[i]>Vp)*(t[i]%%365==0)*( WGAMw[i-1] * Q_W[i-1] / WWw_sum[i-1] * te)
      C_GAM[i] = Q_GAM[i] / WGAMw[i-1]
      Q_W[i] = Q_W[i-1] + (Q_Wcum[i] - Q_Wcum[i-1]) - (V[i]>Vp)*(t[i]%%365==0)*( WGAMw[i-1] * Q_W[i-1] / WWw_sum[i-1] * te)
    }
    
    C_W = Q_W / WWw_sum
    
    ###########################################
    ## EXPORTATION OF RESULTS
    
    
    list(Lw=Lw, WVw=WVw, WVd=WVd, WEw=WEw, WEd=WEd,
         WGAMd=WGAMd, WGAMw=WGAMw, WWd=WWd, WWw_sum=WWw_sum, 
         K_dry=K_dry, dry_content=dry_content,
         WWw_Fulton=WWw_Fulton, WGAMd_prior_repro=WGAMd_prior_repro, WGAMw_prior_repro=WGAMw_prior_repro,
         GSI_WWw_sum=GSI_WWw_sum,GSI_WWw_Fulton=GSI_WWw_Fulton,
         Kw_WWw_sum=Kw_WWw_sum, Kw_WWw_Fulton=Kw_WWw_Fulton,EW_per_gd=EW_per_gd,
         EW_per_gw_WWw_sum=EW_per_gw_WWw_sum, EW_per_gw_WWw_Fulton=EW_per_gw_WWw_Fulton,
         humidity_sum=humidity_sum, humidity_Fulton=humidity_Fulton, rho_V_dry_value=rho_V_dry_value,
         EonWd=EonWd, EonWw_sum=EonWw_sum, EonWw_Fulton=EonWw_Fulton,
         Negg_buffer=Negg_buffer, Negg_buffer_energy=Negg_buffer_energy, Negg_sum=Negg_sum,
         Negg_Fulton=Negg_Fulton,
         C_W=C_W, Q_W=Q_W, Q_GAM=Q_GAM, C_GAM=C_GAM)
    
  }
  
  obs_m <- obs(E=state_m$E, V=state_m$V, Q_Wcum=state_m$Q_Wcum, GAM=state_m$GAM, GAM_prior_repro=state_m$GAM_prior_repro, 
               kappa=kappa_m,Vp=Vp_m,rho_G_dry=rho_Gm_dry, C_G_w=C_Gm_w)
  
  obs_f <- obs(E=state_f$E, V=state_f$V, Q_Wcum=state_f$Q_Wcum, GAM=state_f$GAM, GAM_prior_repro=state_f$GAM_prior_repro,
               kappa=kappa_f,Vp=Vp_f,rho_G_dry=rho_Gf_dry, C_G_w=C_Gf_w)
  
  
  #######################################
  #######   TABLE OF RESULTS   ##########
  #######################################
  
  table <- rbind(
    apply(MARGIN=1,rbind(obs_m$Kw_WWw_sum, obs_m$Kw_WWw_Fulton, obs_m$K_dry,
                         obs_m$EW_per_gw_WWw_sum, obs_m$EW_per_gw_WWw_Fulton, obs_m$EW_per_gd,
                         obs_m$rho_V_dry_value, obs_m$humidity_sum, obs_m$humidity_Fulton,
                         obs_m$EonWd, obs_m$EonWw_sum, obs_m$EonWw_Fulton),
          FUN=mean, na.rm=T),
    apply(MARGIN=1,rbind(obs_f$Kw_WWw_sum, obs_f$Kw_WWw_Fulton, obs_f$K_dry,
                         obs_f$EW_per_gw_WWw_sum, obs_f$EW_per_gw_WWw_Fulton, obs_f$EW_per_gd,
                         obs_f$rho_V_dry_value, obs_f$humidity_sum, obs_f$humidity_Fulton,
                         obs_f$EonWd, obs_f$EonWw_sum, obs_f$EonWw_Fulton),
          FUN=mean, na.rm=T)
  )
  rownames(table) <- c("male","female")
  colnames(table) <- c("Kw_WWw_sum", "Kw_WWw_Fulton", "Kdry", "EW_per_gw_WWw_sum", "EW_per_gw_WWw_Fulton", "EW_per_gd",
                       "rho_V_dry_value", "humidity_sum", "humidity_Fulton", "EonWd", "EonWw_sum", "EonWw_Fulton")
  
  
  #################################
  ##  GRAPHS FOR BIOACCUMULATION ##
  #################################
  
  plot(t,obs_f$Q_GAM,type="l")
  plot(t,obs_f$C_GAM,type="l")
  
# Experiment contamination 3 months / non contamination 3 months (no sex or individual identification)
  if(exp=="1PCB"){
    
    # Physical length
      plot(x=t, y=obs_m$Lw, type="l",col="blue",
           xlab="Time (days)",ylab="Physical length (cm)",
           ylim=c(min(obs_m$Lw,obs_f$Lw,na.rm=T,
                      exp_deconta3_PCB_biometry$"total_length_cm",
                      exp_deconta3_CS_biometry$"total_length_cm"),
                  max(obs_m$Lw,obs_f$Lw,na.rm=T,
                      exp_deconta3_PCB_biometry$"total_length_cm",
                      exp_deconta3_CS_biometry$"total_length_cm")))
      title("Time evolution of physical length",line=+2.5)
      mtext(paste("expérience contamination (3 mois) et décontamination (3 mois)", 
                  ifelse(ingestion==1,"avec nourriture mesurée", "avec f constant")),
            cex=0.8, line=+0.5)
      points(x=t, y=obs_f$Lw, col="red",type="l")
      points(x=n_d1 + exp_deconta3_CS_biometry$"time_d", y=exp_deconta3_CS_biometry$"total_length_cm",pch=3)
      points(x=n_d1 + exp_deconta3_PCB_biometry$"time_d", y=exp_deconta3_PCB_biometry$"total_length_cm",col="orange")
      legend("topleft", lty=c(1,1,0,0), col=c("red","blue","orange","black"), pch=c(NA,NA,1,3),
             legend=c("model simulation for female","model simulation for male",
                      "experimental data for contaminated fish", "experimental data for solvant control fish"))
    
    # Total wet weight
      plot(x=t, y=obs_m$WWw_sum, type="l",col="blue",
           xlab="Time (days)",ylab="Total wet weight (g)",
           ylim=c(min(obs_m$WWw_sum,obs_f$WWw_sum,na.rm=T,
                      exp_deconta3_PCB_biometry$"total_wet_weight_g",
                      exp_deconta3_CS_biometry$"total_wet_weight_g"),
                  max(obs_m$WWw_sum,obs_f$WWw_sum,na.rm=T,
                      exp_deconta3_PCB_biometry$"total_wet_weight_g",
                      exp_deconta3_CS_biometry$"total_wet_weight_g")))
      title("Time evolution of total wet weight",line=+2.5)
      mtext(paste("expérience contamination (3 mois) et décontamination (3 mois)", 
                  ifelse(ingestion==1,"avec nourriture mesurée", "avec f constant")),
            cex=0.8, line=+0.5)
      points(x=t, y=obs_f$WWw_sum, col="red",type="l")
      points(x=n_d1 + exp_deconta3_CS_biometry$"time_d", y=exp_deconta3_CS_biometry$"total_wet_weight_g",pch=3)
    points(x=n_d1 + exp_deconta3_PCB_biometry$"time_d", y=exp_deconta3_PCB_biometry$"total_wet_weight_g",col="orange")
      legend("topleft", lty=c(1,1,0,0), col=c("red","blue","orange","black"), pch=c(NA,NA,1,3),
             legend=c("model simulation for female","model simulation for male",
                      "experimental data for contaminated fish", "experimental data for solvant control fish"))
      
    # Total quantities of contaminant
      plot(t,obs_f$Q_W,type="l",col="red",
           ylab=paste("ng",contaminant), xlab="jours depuis la naissance",
           ylim=c(min(exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant]*as.numeric(as.character(exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)=="total_wet_weight_mean_g"])),
                      obs_f$Q_W,obs_m$Q_W, na.rm=T),
                  max(exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant]*as.numeric(as.character(exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)=="total_wet_weight_mean_g"])),
                      obs_f$Q_W,obs_m$Q_W, na.rm=T)))
      points(t,obs_m$Q_W,type="l",col="blue")
      title(paste("Quantité totale de", contaminant),line=+2.5)
      mtext(paste("\n expérience contamination (3 mois à", C_X,"ng/g de nourriture)",
                  "\n et décontamination (3 mois à", C_X_res,"ng/g de nourriture)", 
                  ifelse(ingestion==1,"avec nourriture mesurée", "avec f constant")),
            cex=0.8, line=+0.5)
      points(x=exp_deconta3_PCB_contam$time_d + n_d1,
             y=exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant]*as.numeric(as.character(exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)=="total_wet_weight_mean_g"])))
      legend("topleft", lty=c(1,0), pch=c(NA,1), legend=c("model simulation","experimental data")) 
      legend("bottomright", pch=19, col=c("red","blue"),legend=c("female","male"))
    
    # Total concentration in contaminant
      plot(t,obs_f$C_W,type="l",col="red",
           ylab=paste("ng",contaminant,"/g frais de poisson"), xlab="jours depuis la naissance",
           ylim=c(min(obs_f$C_W,obs_m$C_W, exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant]),
                  max(obs_f$C_W,obs_m$C_W, exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant])))
      points(t,obs_m$C_W,type="l",col="blue")
      title(paste("Concentration totale en", contaminant),line=+2.5)
      mtext(paste("\n expérience contamination (3 mois à", C_X,"ng/g de nourriture)",
                  "\n et décontamination (3 mois à", C_X_res,"ng/g de nourriture)", 
                  ifelse(ingestion==1,"avec nourriture mesurée", "avec f constant")),
            cex=0.8, line=+0.5)
      points(x=exp_deconta3_PCB_contam$time_d + n_d1,
             y=exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant])
      legend("topleft", lty=c(1,0), pch=c(NA,1), legend=c("model simulation","experimental data")) 
      legend("bottomright", pch=19, col=c("red","blue"),legend=c("female","male"))

    
# Experiment contamination 12 months / non contamination 6 months
      }else {if(exp=="PCBelim"){
        
      # Physical length

        
      # Total wet weight
        
      # Total quantities of contaminant
        plot(t,obs_f$Q_W,type="l",
             ylab=paste("ng",contaminant), xlab="jours depuis la naissance",
             ylim=c(min(exp2_PCB_elim[,colnames(exp2_PCB_elim)==contaminant]*as.numeric(as.character(exp2_PCB_elim[,colnames(exp2_PCB_elim)=="total_wet_weight_g"])),
                        obs_f$Q_W, na.rm=T),
                    max(exp2_PCB_elim[,colnames(exp2_PCB_elim)==contaminant]*as.numeric(as.character(exp2_PCB_elim[,colnames(exp2_PCB_elim)=="total_wet_weight_g"])),
                        obs_f$Q_W, na.rm=T)))
        title(paste("Quantité totale de", contaminant),line=+2.5)
        mtext(paste("\n expérience contamination (12 mois à", C_X,"ng/g de nourriture)",
                    "\n et décontamination (6 mois à", C_X_res,"ng/g de nourriture)", 
                    ifelse(ingestion==1,"avec nourriture mesurée", "avec f constant")),
              cex=0.8, line=+0.5)
        points(x=exp2_PCB_elim$time_d + n_d1,
               y=exp2_PCB_elim[,colnames(exp2_PCB_elim)==contaminant]*as.numeric(as.character(exp2_PCB_elim[,colnames(exp2_PCB_elim)=="total_wet_weight_g"])))
        legend("topleft", lty=c(1,0), pch=c(NA,1), legend=c("model simulation","experimental data")) 
        
      # Total concentration in contaminant
        plot(t,obs_f$C_W,type="l",
             ylab=paste("ng",contaminant,"/g frais de poisson"), xlab="jours depuis la naissance",
             ylim=c(min(obs_f$C_W, exp2_PCB_elim[,colnames(exp2_PCB_elim)==contaminant]),
                    max(obs_f$C_W, exp2_PCB_elim[,colnames(exp2_PCB_elim)==contaminant])))
        title(paste("Concentration totale en", contaminant),line=+2.5)
        mtext(paste("\n expérience contamination (12 mois à", C_X,"ng/g de nourriture)",
                    "\n et décontamination (6 mois à", C_X_res,"ng/g de nourriture)", 
                    ifelse(ingestion==1,"avec nourriture mesurée", "avec f constant")),
              cex=0.8, line=+0.5)
        points(x=exp2_PCB_elim$time_d + n_d1,
               y=exp2_PCB_elim[,colnames(exp2_PCB_elim)==contaminant])
        legend("topleft", lty=c(1,0), pch=c(NA,1), legend=c("model simulation","experimental data")) 
        
        
# Experiment contamination 18 months
      }else {if(exp=="PCBconta18"){
        
      # Physical length
        
      # Total wet weight
        
      # Total quantities of contaminant
        plot(t,obs_f$Q_W,type="l",
             ylab=paste("ng",contaminant), xlab="jours depuis la naissance",
             ylim=c(min(exp2_PCB_conta[,colnames(exp2_PCB_conta)==contaminant]*as.numeric(as.character(exp2_PCB_conta[,colnames(exp2_PCB_conta)=="total_wet_weight_g"])),
                        obs_f$Q_W, na.rm=T),
                    max(exp2_PCB_conta[,colnames(exp2_PCB_conta)==contaminant]*as.numeric(as.character(exp2_PCB_conta[,colnames(exp2_PCB_conta)=="total_wet_weight_g"])),
                        obs_f$Q_W, na.rm=T)))
        title(paste("Quantité totale de", contaminant),line=+2.5)
        mtext(paste("\n expérience contamination (18 mois à", C_X,"ng/g de nourriture) \n",
                    ifelse(ingestion==1,"avec nourriture mesurée", "avec f constant")),
              cex=0.8, line=+0.5)
        points(x=exp2_PCB_conta$time_d + n_d1,
               y=exp2_PCB_conta[,colnames(exp2_PCB_conta)==contaminant]*as.numeric(as.character(exp2_PCB_conta[,colnames(exp2_PCB_conta)=="total_wet_weight_g"])))
        legend("topleft", lty=c(1,0), pch=c(NA,1), legend=c("model simulation","experimental data")) 
      
      # Total concentration in contaminant
        plot(t,obs_f$C_W,type="l",
             ylab=paste("ng",contaminant,"/g frais de poisson"), xlab="jours depuis la naissance",
             ylim=c(min(obs_f$C_W, exp2_PCB_conta[,colnames(exp2_PCB_conta)==contaminant]),
                    max(obs_f$C_W, exp2_PCB_conta[,colnames(exp2_PCB_conta)==contaminant])))
        title(paste("Concentration totale en", contaminant),line=+2.5)
        mtext(paste("\n expérience contamination (18 mois à", C_X,"ng/g de nourriture) \n",
                    ifelse(ingestion==1,"avec nourriture mesurée", "avec f constant")),
              cex=0.8, line=+0.5)
        points(x=exp2_PCB_conta$time_d + n_d1,
               y=exp2_PCB_conta[,colnames(exp2_PCB_conta)==contaminant])
        legend("topleft", lty=c(1,0), pch=c(NA,1), legend=c("model simulation","experimental data")) 
        
        
# Experiment contamination 6 years
      }else {if(exp=="PCBconta6y"){
        
      # Physical length
        
      # Total wet weight
        
      # Total quantities of contaminant
        plot(t,obs_f$Q_W,type="l",col="red",
             ylab=paste("ng",contaminant), xlab="jours depuis la naissance",
             ylim=c(min(exp2_PCB_conta6y[,colnames(exp2_PCB_conta6y)==contaminant]*as.numeric(as.character(exp2_PCB_conta6y[,colnames(exp2_PCB_conta6y)=="total_wet_weight_g"])),
                        obs_f$Q_W, obs_m$Q_W, na.rm=T),
                    max(exp2_PCB_conta6y[,colnames(exp2_PCB_conta6y)==contaminant]*as.numeric(as.character(exp2_PCB_conta6y[,colnames(exp2_PCB_conta6y)=="total_wet_weight_g"])),
                        obs_f$Q_W, obs_m$Q_W, na.rm=T)))
        points(t,obs_m$Q_W,type="l",col="blue")
        title(paste("Quantité totale de", contaminant),line=+2.5)
        mtext(paste("\n expérience contamination (6 ans à", C_X,"ng/g de nourriture) \n",
                    ifelse(ingestion==1,"avec nourriture mesurée", "avec f constant")),
              cex=0.8, line=+0.5)
        points(x=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="I",]$time_d + n_d1,
               y=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="I",colnames(exp2_PCB_conta6y)==contaminant]*as.numeric(as.character(exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="I",colnames(exp2_PCB_conta6y)=="total_wet_weight_g"])),
               col="orange")
        points(x=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="F",]$time_d + n_d1,
               y=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="F",colnames(exp2_PCB_conta6y)==contaminant]*as.numeric(as.character(exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="F",colnames(exp2_PCB_conta6y)=="total_wet_weight_g"])),
               col="red")
        points(x=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="M",]$time_d + n_d1,
               y=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="M",colnames(exp2_PCB_conta6y)==contaminant]*as.numeric(as.character(exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="M",colnames(exp2_PCB_conta6y)=="total_wet_weight_g"])),
               col="blue")
        legend("topleft", lty=c(1,0), pch=c(NA,1), legend=c("model simulation","experimental data"))
        legend("bottomright", pch=19, col=c("red","blue","orange"),
               legend=c("female","male","?"))
        
      # Total concentration in contaminant
        plot(t,obs_f$C_W,type="l",col="red",
             ylab=paste("ng",contaminant,"/g frais de poisson"), xlab="jours depuis la naissance",
             ylim=c(min(obs_f$C_W, exp2_PCB_conta6y[,colnames(exp2_PCB_conta6y)==contaminant],
                        obs_f$C_W,obs_m$C_W,na.rm=T),
                    max(obs_f$C_W, exp2_PCB_conta6y[,colnames(exp2_PCB_conta6y)==contaminant],
                        obs_f$C_W,obs_m$C_W,na.rm=T)))
        points(t,obs_m$C_W,type="l",col="blue")
        title(paste("Concentration totale en", contaminant),line=+2.5)
        mtext(paste("\n expérience contamination (6 ans à", C_X,"ng/g de nourriture) \n",
                    ifelse(ingestion==1,"avec nourriture mesurée", "avec f constant")),
              cex=0.8, line=+0.5)
        points(x=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="I",]$time_d + n_d1,
               y=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="I",colnames(exp2_PCB_conta6y)==contaminant],
               col="orange")
        points(x=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="F",]$time_d + n_d1,
               y=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="F",colnames(exp2_PCB_conta6y)==contaminant],
               col="red")
        points(x=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="M",]$time_d + n_d1,
               y=exp2_PCB_conta6y[exp2_PCB_conta6y$sex=="M",colnames(exp2_PCB_conta6y)==contaminant],
               col="blue")
        legend("topleft", lty=c(1,0), pch=c(NA,1), legend=c("model simulation","experimental data"))
        legend("bottomright", pch=19, col=c("red","blue","orange"),
               legend=c("female","male","?"))
      }else{print("stop")} }}}


  output <- list(t=t, parameters=parameters, state_m=state_m, state_f=state_f, obs_m=obs_m, obs_f=obs_f, table=table)
  
}
  

########################
#    EXP 1 SOLEBEMOL   #
########################

# PCB food=f
CB153_f_1 <- BIOACC(exp="1PCB",contaminant="CB153", C_X=888,log_Kow=6.92,n_d1=6*30,n_d2=84,n_d3=84,rho_food_wet=21612,
                       T_C=19,T_N=NA,ae = 0.8,C_X_res=0,J_food=NA,
                       food_f=0.74, food_m=0.74, ingestion=0)

CB149_f_1 <- BIOACC(exp="1PCB",contaminant="CB149", C_X=420,log_Kow=6.67,n_d1=6*30,n_d2=84,n_d3=84,rho_food_wet=21612,
                T_C=19,T_N=NA,ae = 0.8,C_X_res=0,J_food=NA,
                food_f=0.74, food_m=0.74, ingestion=0)

CB118_f_1 <- BIOACC(exp="1PCB",contaminant="CB118", C_X=454,log_Kow=6.74,n_d1=6*30,n_d2=84,n_d3=84,rho_food_wet=21612,
                T_C=19,T_N=NA,ae = 0.8,C_X_res=0,J_food=NA,
                food_f=0.74, food_m=0.74, ingestion=0)

CB105_f_1 <- BIOACC(exp="1PCB",contaminant="CB105", C_X=228,log_Kow=6.65,n_d1=6*30,n_d2=84,n_d3=84,rho_food_wet=21612,
                T_C=19,T_N=NA,ae = 0.8,C_X_res=0,J_food=NA,
                food_f=0.74, food_m=0.74, ingestion=0)

# PCB food=joules
CB153_J_1 <- BIOACC(exp="1PCB",contaminant="CB153", C_X=888,log_Kow=6.92,n_d1=6*30,n_d2=84,n_d3=84,rho_food_wet=21612,
                    T_C=19,T_N=NA,ae = 0.8,C_X_res=0,J_food=deconta3_PCB_food,
                    food_f=NA, food_m=NA, ingestion=1)

CB149_J_1 <- BIOACC(exp="1PCB",contaminant="CB149", C_X=420,log_Kow=6.67,n_d1=6*30,n_d2=84,n_d3=84,rho_food_wet=21612,
                    T_C=19,T_N=NA,ae = 0.8,C_X_res=0,J_food=deconta3_PCB_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB118_J_1 <- BIOACC(exp="1PCB",contaminant="CB118", C_X=454,log_Kow=6.74,n_d1=6*30,n_d2=84,n_d3=84,rho_food_wet=21612,
                    T_C=19,T_N=NA,ae = 0.8,C_X_res=0,J_food=deconta3_PCB_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB105_J_1 <- BIOACC(exp="1PCB",contaminant="CB105", C_X=228,log_Kow=6.65,n_d1=6*30,n_d2=84,n_d3=84,rho_food_wet=21612,
                    T_C=19,T_N=NA,ae = 0.8,C_X_res=0,J_food=deconta3_PCB_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

PCB1 <- list("CB153_f_1"=CB153_f_1,"CB149_f_1"=CB149_f_1,"CB118_f_1"=CB118_f_1,"CB105_f_1"=CB105_f_1,
             "CB153_J_1"=CB153_J_1,"CB149_J_1"=CB149_J_1,"CB118_J_1"=CB118_J_1,"CB105_J_1"=CB105_J_1)


plot_1PCB <- function(contaminant,n_d1){
  # Physical length
  plot(x=PCB1[[paste(contaminant,"_J_1",sep="")]]$t,
       y=PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$Lw, type="l",col="blue",
       xlab="Time (days)",ylab="Physical length (cm)",
       ylim=c(min(PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$Lw,
                  PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$Lw,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$Lw,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$Lw,
                  na.rm=T,
                  exp_deconta3_PCB_biometry$"total_length_cm",
                  exp_deconta3_CS_biometry$"total_length_cm"),
              max(PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$Lw,
                  PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$Lw,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$Lw,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$Lw,
                  na.rm=T,
                  exp_deconta3_PCB_biometry$"total_length_cm",
                  exp_deconta3_CS_biometry$"total_length_cm")))
  title("Time evolution of physical length",line=+2)
  mtext("experiment of a 3-month contamination period and a 3-month non contamination period",cex=0.8, line=+0.5)
  points(x=PCB1[[paste(contaminant,"_J_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$Lw, col="red",type="l")
  points(x=PCB1[[paste(contaminant,"_f_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$Lw, col="blue",type="l",lty=2)
  points(x=PCB1[[paste(contaminant,"_f_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$Lw, col="red",type="l",lty=2)
  points(x=n_d1 + exp_deconta3_CS_biometry$"time_d", y=exp_deconta3_CS_biometry$"total_length_cm",pch=3)
  points(x=n_d1 + exp_deconta3_PCB_biometry$"time_d", y=exp_deconta3_PCB_biometry$"total_length_cm",col="orange")
  legend("bottomright", lty=c(1,2,0,0,0,0), 
         col=c("black","black","orange","black","red","blue"),
         pch=c(NA,NA,1,3,19,19),
         cex=0.8,
         legend=c("model simulation with measured food ingestion (pX)",
                  "model simulation with fixed food density (f)",
                  "experimental data for contaminated fish", "experimental data for solvant control fish",
                  "female","male"))
  
  # Total wet weight
  plot(x=PCB1[[paste(contaminant,"_J_1",sep="")]]$t,
       y=PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$WWw_sum, type="l",col="blue",
       xlab="Time (days)",ylab="Total wet weight (g)",
       ylim=c(min(PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$WWw_sum,
                  PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$WWw_sum,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$WWw_sum,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$WWw_sum,
                  na.rm=T,
                  exp_deconta3_PCB_biometry$"total_wet_weight_g",
                  exp_deconta3_CS_biometry$"total_wet_weight_g"),
              max(PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$WWw_sum,
                  PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$WWw_sum,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$WWw_sum,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$WWw_sum,
                  na.rm=T,
                  exp_deconta3_PCB_biometry$"total_wet_weight_g",
                  exp_deconta3_CS_biometry$"total_wet_weight_g")))
  title("Time evolution of fish total wet weight",line=+2)
  mtext("experiment of a 3-month contamination period and a 3-month non contamination period",cex=0.8, line=+0.5)
  points(x=PCB1[[paste(contaminant,"_J_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$WWw_sum, col="red",type="l")
  points(x=PCB1[[paste(contaminant,"_f_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$WWw_sum, col="blue",type="l",lty=2)
  points(x=PCB1[[paste(contaminant,"_f_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$WWw_sum, col="red",type="l",lty=2)
  points(x=n_d1 + exp_deconta3_CS_biometry$"time_d", y=exp_deconta3_CS_biometry$"total_wet_weight_g",pch=3)
  points(x=n_d1 + exp_deconta3_PCB_biometry$"time_d", y=exp_deconta3_PCB_biometry$"total_wet_weight_g",col="orange")
  legend("topleft", lty=c(1,2,0,0,0,0), 
         col=c("black","black","orange","black","red","blue"),
         pch=c(NA,NA,1,3,19,19),
         cex=0.8,
         legend=c("model simulation with measured food ingestion (pX)",
                  "model simulation with fixed food density (f)",
                  "experimental data for contaminated fish", "experimental data for solvant control fish",
                  "female","male"))
  
  # Total quantities of contaminant
  plot(x=PCB1[[paste(contaminant,"_J_1",sep="")]]$t,
       y=PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$Q_W, type="l",col="blue",
       ylab=paste("ng",contaminant), xlab="Days from fecundation",
       ylim=c(min(PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$Q_W,
                  PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$Q_W,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$Q_W,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$Q_W,
                  na.rm=T,
                  exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant]*as.numeric(as.character(exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)=="total_wet_weight_mean_g"]))),
              max(PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$Q_W,
                  PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$Q_W,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$Q_W,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$Q_W,
                  na.rm=T,
                  exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant]*as.numeric(as.character(exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)=="total_wet_weight_mean_g"])))))
  title(paste("Total mass of", contaminant, "in fish"),line=+2)
  mtext("experiment of a 3-month contamination period and a 3-month non contamination period",cex=0.8, line=+0.5)
  points(x=PCB1[[paste(contaminant,"_J_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$Q_W, col="red",type="l")
  points(x=PCB1[[paste(contaminant,"_f_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$Q_W, col="blue",type="l",lty=2)
  points(x=PCB1[[paste(contaminant,"_f_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$Q_W, col="red",type="l",lty=2)
  points(x=exp_deconta3_PCB_contam$time_d + n_d1,
         y=exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant]*as.numeric(as.character(exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)=="total_wet_weight_mean_g"])))
  legend("topleft", lty=c(1,2,0,0,0), 
         col=c("black","black","black","red","blue"),
         pch=c(NA,NA,1,19,19),
         cex=0.8,
         legend=c("model simulation with measured food ingestion (pX)",
                  "model simulation with fixed food density (f)",
                  "experimental data","female","male"))
  
  # Total concentration in contaminant
  plot(x=PCB1[[paste(contaminant,"_J_1",sep="")]]$t,
       y=PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$C_W, type="l",col="blue",
       ylab=paste("ng",contaminant,"/gw of fish"), xlab="Days from fecundation",
       ylim=c(min(PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$C_W,
                  PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$C_W,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$C_W,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$C_W,
                  na.rm=T,
                  exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant]),
              max(PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$C_W,
                  PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$C_W,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$C_W,
                  PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$C_W,
                  na.rm=T,
                  exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant])))
  title(paste("Total concentration of", contaminant, "in fish"),line=+2)
  mtext("experiment of 3-month contamination period and 3-month non contamination period",cex=0.8, line=+0.5)
  points(x=PCB1[[paste(contaminant,"_J_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$C_W, col="red",type="l")
  points(x=PCB1[[paste(contaminant,"_f_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$C_W, col="blue",type="l",lty=2)
  points(x=PCB1[[paste(contaminant,"_f_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$C_W, col="red",type="l",lty=2)
  points(x=exp_deconta3_PCB_contam$time_d + n_d1,
         y=exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant])
  legend("topleft", lty=c(1,2,0,0,0), 
         col=c("black","black","black","red","blue"),
         pch=c(NA,NA,1,19,19),
         cex=0.8,
         legend=c("model simulation with measured food ingestion (pX)",
                  "model simulation with fixed food density (f)",
                  "experimental data","female","male"))
}

plot_1PCB("CB153",n_d1=6*30)
plot_1PCB("CB149",n_d1=6*30)
plot_1PCB("CB118",n_d1=6*30)
plot_1PCB("CB105",n_d1=6*30)


########################
#    EXP 2 SOLEBEMOL   #
########################


# PCB ELIM food=f
CB153_f_2 <- BIOACC(exp="PCBelim",contaminant="CB153", C_X_res=0, C_X=300,log_Kow=6.92,
                    n_d1=6*30,n_d2=12*30,n_d3=6*30,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

CB149_f_2 <- BIOACC(exp="PCBelim",contaminant="CB149", C_X_res=0, C_X=150,log_Kow=6.67,
                    n_d1=6*30,n_d2=12*30,n_d3=6*30,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

CB118_f_2 <- BIOACC(exp="PCBelim",contaminant="CB118", C_X_res=0, C_X=150,log_Kow=6.74,
                    n_d1=6*30,n_d2=12*30,n_d3=6*30,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

CB105_f_2 <- BIOACC(exp="PCBelim",contaminant="CB105", C_X_res=0, C_X=75,log_Kow=6.65,
                    n_d1=6*30,n_d2=12*30,n_d3=6*30,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

# PCB ELIM food=joules
CB153_J_2 <- BIOACC(exp="PCBelim",contaminant="CB153", C_X_res=0, C_X=300,log_Kow=6.92,
                    n_d1=6*30,n_d2=12*30,n_d3=6*30,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB149_J_2 <- BIOACC(exp="PCBelim",contaminant="CB149", C_X_res=0, C_X=150,log_Kow=6.67,
                    n_d1=6*30,n_d2=12*30,n_d3=6*30,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB118_J_2 <- BIOACC(exp="PCBelim",contaminant="CB118", C_X_res=0, C_X=150,log_Kow=6.74,
                    n_d1=6*30,n_d2=12*30,n_d3=6*30,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB105_J_2 <- BIOACC(exp="PCBelim",contaminant="CB105", C_X_res=0, C_X=75,log_Kow=6.65,
                    n_d1=6*30,n_d2=12*30,n_d3=6*30,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

# PCB conta18 food=f
CB153_f_3 <- BIOACC(exp="PCBconta18",contaminant="CB153", C_X_res=0, C_X=300,log_Kow=6.92,
                    n_d1=6*30,n_d2=18*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

CB149_f_3 <- BIOACC(exp="PCBconta18",contaminant="CB149", C_X_res=0, C_X=150,log_Kow=6.67,
                    n_d1=6*30,n_d2=18*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

CB118_f_3 <- BIOACC(exp="PCBconta18",contaminant="CB118", C_X_res=0, C_X=150,log_Kow=6.74,
                    n_d1=6*30,n_d2=18*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

CB105_f_3 <- BIOACC(exp="PCBconta18",contaminant="CB105", C_X_res=0, C_X=75,log_Kow=6.65,
                    n_d1=6*30,n_d2=18*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

# PCB conta18 food=joules
CB153_J_3 <- BIOACC(exp="PCBconta18",contaminant="CB153", C_X_res=0, C_X=300,log_Kow=6.92,
                    n_d1=6*30,n_d2=18*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB149_J_3 <- BIOACC(exp="PCBconta18",contaminant="CB149", C_X_res=0, C_X=150,log_Kow=6.67,
                    n_d1=6*30,n_d2=18*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB118_J_3 <- BIOACC(exp="PCBconta18",contaminant="CB118", C_X_res=0, C_X=150,log_Kow=6.74,
                    n_d1=6*30,n_d2=18*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB105_J_3 <- BIOACC(exp="PCBconta18",contaminant="CB105", C_X_res=0, C_X=75,log_Kow=6.65,
                    n_d1=6*30,n_d2=18*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

PCB2 <- list("CB153_f_2"=CB153_f_2,"CB149_f_2"=CB149_f_2,"CB118_f_2"=CB118_f_2,"CB105_f_2"=CB105_f_2,
             "CB153_J_2"=CB153_J_2,"CB149_J_2"=CB149_J_2,"CB118_J_2"=CB118_J_2,"CB105_J_2"=CB105_J_2,
             "CB153_f_3"=CB153_f_3,"CB149_f_3"=CB149_f_3,"CB118_f_3"=CB118_f_3,"CB105_f_3"=CB105_f_3,
             "CB153_J_3"=CB153_J_3,"CB149_J_3"=CB149_J_3,"CB118_J_3"=CB118_J_3,"CB105_J_3"=CB105_J_3)

plot_2PCB <- function(contaminant,n_d1){
# Total quantities of contaminant
plot(x=PCB2[[paste(contaminant,"_J_2",sep="")]]$t,
     y=PCB2[[paste(contaminant,"_J_2",sep="")]][["obs_m"]]$Q_W, type="l",col="blue",
     ylab=paste("ng",contaminant), xlab="Days from fecundation",
     ylim=c(min(PCB2[[paste(contaminant,"_J_2",sep="")]][["obs_m"]]$Q_W,
                PCB2[[paste(contaminant,"_J_2",sep="")]][["obs_f"]]$Q_W,
                PCB2[[paste(contaminant,"_J_3",sep="")]][["obs_m"]]$Q_W,
                PCB2[[paste(contaminant,"_J_3",sep="")]][["obs_f"]]$Q_W,
                na.rm=T,
                exp2_PCB_contamination[,colnames(exp2_PCB_contamination)==contaminant]*as.numeric(as.character(exp2_PCB_contamination[,colnames(exp2_PCB_contamination)=="total_wet_weight_g"]))),
            max(PCB2[[paste(contaminant,"_J_2",sep="")]][["obs_m"]]$Q_W,
                PCB2[[paste(contaminant,"_J_2",sep="")]][["obs_f"]]$Q_W,
                PCB2[[paste(contaminant,"_J_3",sep="")]][["obs_m"]]$Q_W,
                PCB2[[paste(contaminant,"_J_3",sep="")]][["obs_f"]]$Q_W,
                na.rm=T,
                exp2_PCB_contamination[,colnames(exp2_PCB_contamination)==contaminant]*as.numeric(as.character(exp2_PCB_contamination[,colnames(exp2_PCB_contamination)=="total_wet_weight_g"])))))
title(paste("Total mass of", contaminant, "in fish"),line=+2)
mtext("experiment of a 18-month contamination period VS 1 year contamination period and a 6-month non contamination period",cex=0.8, line=+0.5)
points(x=PCB2[[paste(contaminant,"_J_2",sep="")]]$t, y=PCB2[[paste(contaminant,"_J_2",sep="")]][["obs_f"]]$Q_W, col="red",type="l")
points(x=PCB2[[paste(contaminant,"_J_3",sep="")]]$t, y=PCB2[[paste(contaminant,"_J_3",sep="")]][["obs_m"]]$Q_W, col="blue",type="l",lty=2)
points(x=PCB2[[paste(contaminant,"_J_3",sep="")]]$t, y=PCB2[[paste(contaminant,"_J_3",sep="")]][["obs_f"]]$Q_W, col="red",type="l",lty=2)
points(x=exp2_PCB_elim$time_d + n_d1,pch=1,
       y=exp2_PCB_elim[,colnames(exp2_PCB_elim)==contaminant]*as.numeric(as.character(exp2_PCB_elim[,colnames(exp2_PCB_elim)=="total_wet_weight_mean_g"])))
points(x=exp2_PCB_conta$time_d + n_d1,pch=2,
       y=exp2_PCB_conta[,colnames(exp2_PCB_conta)==contaminant]*as.numeric(as.character(exp2_PCB_conta[,colnames(exp2_PCB_conta)=="total_wet_weight_mean_g"])))

legend("topleft", lty=c(1,2,0,0,0), 
       col=c("black","black","black","red","blue"),
       pch=c(NA,NA,1,2,19,19),
       cex=0.8,
       legend=c("model simulation for a 18-months contamination period",
                "model simulation for a 1 year contamination period and a 6-month non contamination period",
                "experimental data elim","experimental data conta","female","male"))
}

plot_2PCB("CB153",n_d1=6*30)
plot_1PCB("CB149",n_d1=6*30)
plot_1PCB("CB118",n_d1=6*30)
plot_1PCB("CB105",n_d1=6*30)


# Total concentration in contaminant
plot(x=PCB1[[paste(contaminant,"_J_1",sep="")]]$t,
     y=PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$C_W, type="l",col="blue",
     ylab=paste("ng",contaminant,"/gw of fish"), xlab="Days from fecundation",
     ylim=c(min(PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$C_W,
                PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$C_W,
                PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$C_W,
                PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$C_W,
                na.rm=T,
                exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant]),
            max(PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_m"]]$C_W,
                PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$C_W,
                PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$C_W,
                PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$C_W,
                na.rm=T,
                exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant])))
title(paste("Total concentration of", contaminant, "in fish"),line=+2)
mtext("experiment of 3-month contamination period and 3-month non contamination period",cex=0.8, line=+0.5)
points(x=PCB1[[paste(contaminant,"_J_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_J_1",sep="")]][["obs_f"]]$C_W, col="red",type="l")
points(x=PCB1[[paste(contaminant,"_f_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_m"]]$C_W, col="blue",type="l",lty=2)
points(x=PCB1[[paste(contaminant,"_f_1",sep="")]]$t, y=PCB1[[paste(contaminant,"_f_1",sep="")]][["obs_f"]]$C_W, col="red",type="l",lty=2)
points(x=exp_deconta3_PCB_contam$time_d + n_d1,
       y=exp_deconta3_PCB_contam[,colnames(exp_deconta3_PCB_contam)==contaminant])
legend("topleft", lty=c(1,2,0,0,0), 
       col=c("black","black","black","red","blue"),
       pch=c(NA,NA,1,19,19),
       cex=0.8,
       legend=c("model simulation with measured food ingestion (pX)",
                "model simulation with fixed food density (f)",
                "experimental data","female","male"))
}


# PCB conta6y food=f
CB153_f_4 <- BIOACC(exp="PCBconta6y",contaminant="CB153", C_X_res=0, C_X=300,log_Kow=6.92,
                    n_d1=6*30,n_d2=6*12*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

CB149_f_4 <- BIOACC(exp="PCBconta6y",contaminant="CB149", C_X_res=0, C_X=150,log_Kow=6.67,
                    n_d1=6*30,n_d2=6*12*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

CB118_f_4 <- BIOACC(exp="PCBconta6y",contaminant="CB118", C_X_res=0, C_X=150,log_Kow=6.74,
                    n_d1=6*30,n_d2=6*12*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

CB105_f_4 <- BIOACC(exp="PCBconta6y",contaminant="CB105", C_X_res=0, C_X=75,log_Kow=6.65,
                    n_d1=6*30,n_d2=6*12*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=NA,
                    food_f=0.74, food_m=0.74, ingestion=0)

# PCB conta6y food=joules
CB153_J_4 <- BIOACC(exp="PCBconta6y",contaminant="CB153", C_X_res=0, C_X=300,log_Kow=6.92,
                    n_d1=6*30,n_d2=6*12*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB149_J_4 <- BIOACC(exp="PCBconta6y",contaminant="CB149", C_X_res=0, C_X=150,log_Kow=6.67,
                    n_d1=6*30,n_d2=6*12*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB118_J_4 <- BIOACC(exp="PCBconta6y",contaminant="CB118", C_X_res=0, C_X=150,log_Kow=6.74,
                    n_d1=6*30,n_d2=6*12*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)

CB105_J_4 <- BIOACC(exp="PCBconta6y",contaminant="CB105", C_X_res=0, C_X=75,log_Kow=6.65,
                    n_d1=6*30,n_d2=6*12*30,n_d3=0,rho_food_wet=21612,
                    T_C=0,T_N=exp2_temperature,ae = 0.8,J_food=exp2_food,
                    food_f=0.74, food_m=0.74, ingestion=1)