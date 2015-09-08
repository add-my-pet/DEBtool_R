## mydata_iso_221_var
# created 2011/05/04 by Bas Kooijman, modified 2012/01/30
# isomorph with 1 structure
# 2 reserves: protein (1) and non-protein (2)
#  somatic maintenance and growth overhead preferably paid from non-protein
# 2 types of food
#  preference depends on stress of non-filled reserve
library(deSolve)
source('DEBtoolR/iso221/iso_21_var.R'); source('DEBtoolR/iso221/iso_21_b_var.R'); source('DEBtoolR/iso221/sgr_iso_21_var.R');
source('DEBtoolR/iso221/iso_221_var.R'); source('DEBtoolR/iso221/diso_221b_var.R')


## set parameters at T_ref = 293 K
M_X1      = 1e-3;   M_X2      = 1e-3;  # mol, size of food particle of type i
F_X1m     = 10;     F_X2m     = 10;    # dm^2/d.cm^2, {F_Xim} spec searching rates
y_P1X1    = 0.15;   y_P2X2    = 0.15;  # mol/mol, yield of feaces i on food i
y_E1X1    = 0.8*0.3;  y_E2X1= 0.8*0.7;  # mol/mol, yield of reserve Ei on food X1 (protein, non-protein)
y_E1X2    = 0.8*0.7;   y_E2X2= 0.8*0.3;  # mol/mol, yield of reserve Ei on food X2 (protein, non-protein)
J_X1Am    = 2.0e-3; J_X2Am    = 2.0e-3;# mol/d.cm^2, {J_XiAm} max specific ingestion rate for food Xi
v         = 0.05;   kap       = 0.8;   # cm/d, energy conductance, 
                                       # -, allocation fraction to soma
mu_E1     = 4e5;    mu_E2     = 6e5;   # J/mol, chemical potential of reserve i
mu_V      = 5e5;    j_E1M     = 0.09;  # J/mol, chemical potenial of structure; j_E2M=j_E1M * mu_E1/ mu_E2
                                       # mol/d.mol, specific som maint costs
J_E1T     = 0;      MV        = 4e-3;  # mol/d.cm^2, {J_E1T}, spec surface-area-linked som maint costs J_E1T/ J_E2T = j_E1M/ j_E2M
                                       # mol/cm^3, [M_V] density of structure
k_J       = 0.002;  k1_J      = 0.002; # 1/d, mat maint rate coeff, spec rejuvenation rate                                    
kap_G      = 0.8;   del_V     = 0.8;   # -, growth efficiency
                                       # -, threshold for death by  shrinking
kap_E1    = 0.99;   kap_E2    = 0.99;     # -, fraction of rejected mobilised flux that is returned to reserve
# since j_E1P = 0, kap_E1 is not relevant
kap_R1    = 0.95;   kap_R2    = 0.95;  # -, reproduction efficiency for reserve i
E_Hb      = 1e1;    E_Hp      = 2e4;   # J, maturity thresholds at birth, puberty
T_A       = 8000;   h_H       = 1e-5;  # K, Arrhenius temperature
                                       # 1/d, hazerd due to rejuvenation
h_a       = 2e-8;   s_G       = 1e-4;  # 1/d^2, aging acceleration
                                       # -, Gompertz stress coefficient

# pack parameters
par_iso_221 = c(
#      1       2       3       4       5       6       7       8
    M_X1,   M_X2,  F_X1m,  F_X2m, y_P1X1, y_P2X2, y_E1X1, y_E2X1,
#      9      10      11      12      13      14      15      16
  y_E1X2, y_E2X2, J_X1Am, J_X2Am,      v,    kap,  mu_E1,  mu_E2,
#     17      18      19      20      21      22      23      24 
    mu_V,  j_E1M,  J_E1T,     MV,    k_J,   k1_J,  kap_G,  del_V, 
#     25      26      27      28      29      30      31      32   
  kap_E1, kap_E2, kap_R1, kap_R2,   E_Hb,   E_Hp,    T_A,    h_H, 
#     33      34 
     h_a,    s_G);

# set chemical indices
#    X1   X2    V   E1   E2   P1   P2  organics
n_O = c(
      1,    1,    1,    1,    1,    1,    1,    # C
      1.8,  1.8,  1.8,  1.61, 2.0,  1.8,  1.8,  # H
      0.5,  0.5,  0.5,  0.33, 0.6,  0.5,  0.6,  # O
      0.2,  0.2,  0.2,  0.28, 0.0,  0.2,  0.0); # N
  
#     C    H    O    N                 minerals
n_M = c(
      1,    0,    0,    0, # C
      0,    2,    0,    3, # H
      2,    1,    2,    0, # O
      0,    0,    0,    1);# N

## set environmental variables
t = seq(0,100,length=11); tXT = cbind(t, t, t, t); # d, time points
tXT[,2] = 20000;     tXT[,3] = 20000;               # mol/dm^2, food densities (don't need to be constant)
tXT[,4] = 293;                               # K, temperature (does not need to be constant)

## get state at birth
iso_21_b_var_out = iso_21_b_var(tXT, par_iso_221);

var_b=unlist(iso_21_b_var_out[1]); a_b=unlist(iso_21_b_var_out[2]); M_E10=unlist(iso_21_b_var_out[3]); M_E20=unlist(iso_21_b_var_out[4]); 

## run iso_221
iso_221_var_out = iso_221_var(tXT, var_b, par_iso_221, n_O, n_M); # from birth to t = tXT(end,1)
var<-matrix(data=unlist(iso_221_var_out[1]), nrow = nrow(tXT), ncol = length(var_b)); flux<-matrix(data=unlist(iso_221_var_out[2]), nrow = nrow(tXT), ncol = 6);

GFout = matrix(data=0,nrow=11,ncol=34);
for(i in 1:11){
GFout[i,] = t(as.matrix(unlist(diso_221b_var(tXT[i,1], var[i,], tXT, par_iso_221))));
}

# GF experiment 1, fixed diets

fract<-c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99)
for(k in 1:11){
# continue with a period with only food type 1
y_E1X1    = 0.8*fract[k];  y_E2X1= 0.8*(1-fract[k]);  # mol/mol, yield of reserve Ei on food X1 (protein, non-protein)
y_E1X2    = 0.8*(1-fract[k]);   y_E2X2= 0.8*fract[k];  # mol/mol, yield of reserve Ei on food X2 (protein, non-protein)

# pack parameters
par_iso_221 = c(
#      1       2       3       4       5       6       7       8
    M_X1,   M_X2,  F_X1m,  F_X2m, y_P1X1, y_P2X2, y_E1X1, y_E2X1,
#      9      10      11      12      13      14      15      16
  y_E1X2, y_E2X2, J_X1Am, J_X2Am,      v,    kap,  mu_E1,  mu_E2,
#     17      18      19      20      21      22      23      24 
    mu_V,  j_E1M,  J_E1T,     MV,    k_J,   k1_J,  kap_G,  del_V, 
#     25      26      27      28      29      30      31      32   
  kap_E1, kap_E2, kap_R1, kap_R2,   E_Hb,   E_Hp,    T_A,    h_H, 
#     33      34 
     h_a,    s_G);


# continue with a period with only food type 2
t2 = seq(100,10e3,length=199); tXT2 = cbind(t2, t2, t2, t2); # d, set time points
tXT2[,2] = 2e4; tXT2[,3] = 0; tXT2[,4] = 293;         # set food, temp
var_0 = var[nrow(var),];                                   # copy last state to initial state
iso_221_var_out  = iso_221_var(tXT2, var_0, par_iso_221, n_O, n_M); # run iso_221_var
var2=matrix(data=unlist(iso_221_var_out[1]), nrow = nrow(tXT2), ncol = length(var_b)); flux2=matrix(data=unlist(iso_221_var_out[2]), nrow = nrow(tXT2), ncol = 6);

GFout2 = matrix(data=0,nrow=199,ncol=34);
for(i in 1:199){
GFout2[i,] = t(as.matrix(unlist(diso_221b_var(tXT2[i,1], var2[i,], tXT2, par_iso_221))));
}

# catenate results for plotting
t3 = c(t, t2); var3 = rbind(var, var2); flux3 = rbind(flux, flux2); GFout3 = rbind(GFout, GFout2); tXT3 = rbind(tXT, tXT2);

## plot results
# unpack var: (n,13)-matrix with variables
#  cM_X1, cM_X2, M_E1, M_E2, M_V, M_H, cM_ER1, cM_ER2, q, h, S
#    cum food eaten, reserves, (max)structure, (max)maturity , cum allocation to reprod, accel, hazard, surv
 cM_X1 = var3[, 1]; cM_X2   = var3[, 2]; # mol, cumulative ingested food
 M_E1  = var3[, 3]; M_E2    = var3[, 4]; # mol, reserve
 E_H   = var3[, 5]; max_E_H = var3[, 6]; # J, maturity, max maturity
 M_V   = var3[, 7]; max_M_V = var3[, 8]; # mol, structure, max structure
 cM_E1R= var3[, 9]; cM_E2R  = var3[,10]; # mol, cumulative reprod
 q     = var3[,11]; h       = var3[,12]; # 1/d^2, 1/d, aging acceleration, hazard
 S     = var3[,13];                      # -, survival probability

# unpack flux: (n,20)-matrix with fluxes (most of it still needs to be coded)
#  f1, f2, J_X1A, J_X2A, J_E1A, J_E2A, J_EC1, J_EC2, J_EM1, J_EM2, J_VG, ...
#  J_E1J, J_E2J, J_E1R, J_E2R, R, ...
#  J_C, J_H, J_O, J_N
#    func responses, food eaten, assim, mobilisation, som. maint, growth, ...
#    mat. maint, maturation, reprod rate, ...
#    CO2, H20, O2, NH3
 f1 = flux3[,1];   f2 = flux3[,2];           # -, scaled functional response
 s1 = flux3[,3];   s2 = flux3[,4];           # -, stress coefficients
 rho_X1X2 = flux3[,5]; rho_X2X1 = flux3[,6]; # -, competition coefficients

drymass=var3[,3]+var3[,4]+var3[,7];
repro=min(var3[210,9],var3[210,10]);
OUT = cbind(tXT3, var3, flux3, drymass, GFout3); 
colnames(OUT) = c('Time', 'X1', 'X2', 'Tb', 'cM_X1', 'cM_X2', 'M_E1', 'M_E2', 'E_H', 'max_E_H', 'M_V', 'max_M_V', 'cM_ER1', 'cM_ER2', 'q', 'h', 'surviv', 'f1', 'f2', 's1', 's2', 'rho_X1X2', 'rho_X2X1', 'drymass', 'pyield','cyield','ingest1_(I1)','ingest2_(I2)','assim1_(I1-De1)','assim2_(I2-De2)','mobil1_(I1-De1-Rs1+Rs1(part_of))','mobil2_(I2-De2-Rs2+Rs2(part_of))','grow_tiss1_(Rg1)','grow_tiss2_(Rg2)','grow_ohead1_(Dmp1a)','grow_ohead2_(Dmp2a)','maint1_(Dmb1_+_Dma1)','maint2_(Dmb2_+_Dma2)','mat_maint1_(Dmp1b)','mat_maint2_(Dmp2b)','maturation1_(Dmp1c)','maturation2_(Dmp2c)','repro_tiss1_(Rr1)','repro_tiss2_(Rr2)','repro_ohead1_(Dmp1d)','repro_ohead2_(Dmp2d)','lost_rejres1_(Dc1)','lost_rejres2_(Dc2)','ret_rejres1_(Rs1(part_of))','ret_rejres2_(Rs2(part_of))','dres1_(Rs1)','dres2_(Rs2)','faeces1_(De1)','faeces2_(De2)','tot_retained1_(R1)','tot_retained2_(R2)','growth_repro1_(Rg+Rr)','growth_repro2__(Rg+Rr)'); 

if(k==1){
  repro_fixed<-repro
  GF_fixed<-as.data.frame(OUT)
}else{
  GF_fixed<-rbind(GF_fixed,as.data.frame(OUT))
  repro_fixed<-c(repro_fixed,repro)
}
  
} # end loop through food types

plot(repro_fixed~fract,ylim=c(0,30))

# GF experiment 2, choice

fract<-c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99)
for(k in 1:11){
# continue with a period with only food type 1
y_E1X1    = 0.8*fract[k];  y_E2X1= 0.8*(1-fract[k]);  # mol/mol, yield of reserve Ei on food X1 (protein, non-protein)
y_E1X2    = 0.8*(1-fract[k]);   y_E2X2= 0.8*fract[k];  # mol/mol, yield of reserve Ei on food X2 (protein, non-protein)

# pack parameters
par_iso_221 = c(
#      1       2       3       4       5       6       7       8
    M_X1,   M_X2,  F_X1m,  F_X2m, y_P1X1, y_P2X2, y_E1X1, y_E2X1,
#      9      10      11      12      13      14      15      16
  y_E1X2, y_E2X2, J_X1Am, J_X2Am,      v,    kap,  mu_E1,  mu_E2,
#     17      18      19      20      21      22      23      24 
    mu_V,  j_E1M,  J_E1T,     MV,    k_J,   k1_J,  kap_G,  del_V, 
#     25      26      27      28      29      30      31      32   
  kap_E1, kap_E2, kap_R1, kap_R2,   E_Hb,   E_Hp,    T_A,    h_H, 
#     33      34 
     h_a,    s_G);


# continue with a period with only food type 2
t2 = seq(100,10e3,length=199); tXT2 = cbind(t2, t2, t2, t2); # d, set time points
tXT2[,2] = 2e4; tXT2[,3] = 2e4; tXT2[,4] = 293;         # set food, temp
var_0 = var[nrow(var),];                                   # copy last state to initial state
iso_221_var_out  = iso_221_var(tXT2, var_0, par_iso_221, n_O, n_M); # run iso_221_var
var2=matrix(data=unlist(iso_221_var_out[1]), nrow = nrow(tXT2), ncol = length(var_b)); flux2=matrix(data=unlist(iso_221_var_out[2]), nrow = nrow(tXT2), ncol = 6);

GFout2 = matrix(data=0,nrow=199,ncol=34);
for(i in 1:199){
GFout2[i,] = t(as.matrix(unlist(diso_221b_var(tXT2[i,1], var2[i,], tXT2, par_iso_221))));
}

# catenate results for plotting
t3 = c(t, t2); var3 = rbind(var, var2); flux3 = rbind(flux, flux2); GFout3 = rbind(GFout, GFout2); tXT3 = rbind(tXT, tXT2);

## plot results
# unpack var: (n,13)-matrix with variables
#  cM_X1, cM_X2, M_E1, M_E2, M_V, M_H, cM_ER1, cM_ER2, q, h, S
#    cum food eaten, reserves, (max)structure, (max)maturity , cum allocation to reprod, accel, hazard, surv
 cM_X1 = var3[, 1]; cM_X2   = var3[, 2]; # mol, cumulative ingested food
 M_E1  = var3[, 3]; M_E2    = var3[, 4]; # mol, reserve
 E_H   = var3[, 5]; max_E_H = var3[, 6]; # J, maturity, max maturity
 M_V   = var3[, 7]; max_M_V = var3[, 8]; # mol, structure, max structure
 cM_E1R= var3[, 9]; cM_E2R  = var3[,10]; # mol, cumulative reprod
 q     = var3[,11]; h       = var3[,12]; # 1/d^2, 1/d, aging acceleration, hazard
 S     = var3[,13];                      # -, survival probability

# unpack flux: (n,20)-matrix with fluxes (most of it still needs to be coded)
#  f1, f2, J_X1A, J_X2A, J_E1A, J_E2A, J_EC1, J_EC2, J_EM1, J_EM2, J_VG, ...
#  J_E1J, J_E2J, J_E1R, J_E2R, R, ...
#  J_C, J_H, J_O, J_N
#    func responses, food eaten, assim, mobilisation, som. maint, growth, ...
#    mat. maint, maturation, reprod rate, ...
#    CO2, H20, O2, NH3
 f1 = flux3[,1];   f2 = flux3[,2];           # -, scaled functional response
 s1 = flux3[,3];   s2 = flux3[,4];           # -, stress coefficients
 rho_X1X2 = flux3[,5]; rho_X2X1 = flux3[,6]; # -, competition coefficients

drymass=var3[,3]+var3[,4]+var3[,7];
repro=min(var3[210,9],var3[210,10]);
OUT = cbind(tXT3, var3, flux3, drymass, GFout3); 
colnames(OUT) = c('Time', 'X1', 'X2', 'Tb', 'cM_X1', 'cM_X2', 'M_E1', 'M_E2', 'E_H', 'max_E_H', 'M_V', 'max_M_V', 'cM_ER1', 'cM_ER2', 'q', 'h', 'surviv', 'f1', 'f2', 's1', 's2', 'rho_X1X2', 'rho_X2X1', 'drymass', 'pyield','cyield','ingest1_(I1)','ingest2_(I2)','assim1_(I1-De1)','assim2_(I2-De2)','mobil1_(I1-De1-Rs1+Rs1(part_of))','mobil2_(I2-De2-Rs2+Rs2(part_of))','grow_tiss1_(Rg1)','grow_tiss2_(Rg2)','grow_ohead1_(Dmp1a)','grow_ohead2_(Dmp2a)','maint1_(Dmb1_+_Dma1)','maint2_(Dmb2_+_Dma2)','mat_maint1_(Dmp1b)','mat_maint2_(Dmp2b)','maturation1_(Dmp1c)','maturation2_(Dmp2c)','repro_tiss1_(Rr1)','repro_tiss2_(Rr2)','repro_ohead1_(Dmp1d)','repro_ohead2_(Dmp2d)','lost_rejres1_(Dc1)','lost_rejres2_(Dc2)','ret_rejres1_(Rs1(part_of))','ret_rejres2_(Rs2(part_of))','dres1_(Rs1)','dres2_(Rs2)','faeces1_(De1)','faeces2_(De2)','tot_retained1_(R1)','tot_retained2_(R2)','growth_repro1_(Rg+Rr)','growth_repro2__(Rg+Rr)'); 

if(k==1){
  repro_choice<-repro
  GF_choice<-as.data.frame(OUT)
}else{
  GF_choice<-rbind(GF_choice,as.data.frame(OUT))
  repro_choice<-c(repro_choice,repro)
}
  
} # end loop through food types

points(repro_choice~fract,col='blue')

# for a reason I don't understand, the subset below by the variable 'fract' fails for some values, but writing
# and reading back in resolves it
write.csv(GF_choice,'GF_choice.csv')
write.csv(GF_fixed,'GF_fixed.csv')
GF_choice<-read.csv('GF_choice.csv')
GF_fixed<-read.csv('GF_fixed.csv')
file.remove('GF_choice.csv')
file.remove('GF_fixed.csv')

GF_choice1=subset(GF_choice,Time>100)
GF_fixed1=subset(GF_fixed,Time>100)

agg_choice_max<-aggregate(GF_choice1, by = list(GF_choice1$pyield),FUN = max)
agg_choice_sum<-aggregate(GF_choice1, by = list(GF_choice1$pyield),FUN = sum)

agg_fixed_max<-aggregate(GF_fixed1, by = list(GF_fixed1$pyield),FUN = max)
agg_fixed_sum<-aggregate(GF_fixed1, by = list(GF_fixed1$pyield),FUN = sum)

for(i in 1:11){
  choice<-subset(GF_choice1,pyield==fract[i])
  if(i==1){
  with(choice,plot(M_V~Time, type='l',col=i,ylim=c(0,max(GF_choice1$M_V))))
  }else{
    with(choice,points(M_V~Time, type='l',col=i,ylim=c(0,max(GF_choice1$M_V))))
}
}

for(i in 1:11){
  fixed<-subset(GF_fixed1,GF_fixed1$pyield==fract[i])
  if(i==1){  
   with(fixed,plot(M_V~Time, type='l',col=i,ylim=c(0,max(GF_fixed1$M_V))))
  }else{
   with(fixed,points(M_V~Time, type='l',col=i,ylim=c(0,max(GF_fixed1$M_V))))
  }
}








# figures
par(mfrow = c(2,4)) # set up for 8 plots in 2 rows
plot(cM_X1~t3,col='blue',type='l',ylim=c(0,max(max(cM_X1),max(cM_X2))), xlab='time since birth, d', ylab='cum food eaten, mol')
points(cM_X2~t3,col='red',type='l')
legend(0,80,c(expression('X'[1]),expression('X'[2])),lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"),bty = "n") 

plot(f1~t3,col='blue',type='l',ylim=c(0,1), xlab='time since birth, d', ylab='scaled func resp, -')
points(f2~t3,col='red',type='l')

plot(M_E1~t3,col='blue',type='l',ylim=c(0,max(max(M_E1),max(M_E2))), xlab='time since birth, d', ylab='reserve, mol')
points(M_E2~t3,col='red',type='l')
legend(2000,0.3,c('protein','non-protein'),lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"),bty = "n") 

plot(M_E1/ M_V~t3,col='blue',type='l',ylim=c(0,max(max(M_E1/ M_V),max(M_E2/ M_V))), xlab='time since birth, d', ylab='reserve density, mol/mol')
points(M_E2/ M_V~t3,col='red',type='l')

plot(M_V~t3, col='green',type='l', xlab='time since birth, d', ylab='structure, mol')

plot((M_V/ MV)^(1/3)~t3, col='green',type='l', xlab='time since birth, d', ylab='length, cm')

plot(E_H~t3, col='green',type='l', xlab='time since birth, d', ylab='maturity, J')

plot(cM_E1R~t3, col='blue',type='l', ylim=c(0,max(max(cM_E1R),max(cM_E1R))), xlab='time since birth, d', ylab='cum reprod, mol')
points(cM_E2R~t3, col='red',type='l')


par(mfrow = c(1,3)) # set up for 3 plots in 1 rows
plot(s1~t3,col='blue',type='l',ylim=c(0,max(max(s1),max(s2))), xlab='time since birth, d', ylab='stress coefficients, -')
points(s2~t3,col='red',type='l')
legend(2000,0.1,c('protein','non-protein'),lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"),bty = "n") 

plot(rho_X1X2~t3,col='blue',type='l',ylim=c(0,max(max(rho_X1X2),max(rho_X2X1))), xlab='time since birth, d', ylab='competition coefficients, -')
points(rho_X2X1~t3,col='red',type='l')
legend(0,3,c(expression(rho['X'[1]*'X'[2]]),expression(rho['X'[1]*'X'[2]])),lty=c(1,1), lwd=c(2.5,2.5),col=c("blue","red"),bty = "n") 

plot(S~t3, col='green',type='l', xlab='time since birth, d', ylab='survival prob, -')

