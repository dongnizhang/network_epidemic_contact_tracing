
library(metR)
#install.packages('metR')
library(ggplot2)
#install.packages("pracma")
library(pracma)#AUC
library(rootSolve)
#install.packages("sparsevar")
library(sparsevar)
library(cowplot)
require(gridExtra)

######################################################
## Basic reproduction number: without CT
## assume the degree D follows a poisson distribution.
######################################################

R_0 = function(mu, beta_N, beta_G, tau){
  
# For type-2 (infected through global contacts), the expected size-biased degree minus 1 would be mu
D_tilde_minus_1 = mu

############################################  
# elements of the next generation matrix
############################################
mNN = D_tilde_minus_1*(1-exp(-beta_N*tau))
mHN = mu*(1-exp(-beta_N*tau))
mHH = beta_G*tau
mNH = beta_G*tau

M <- matrix(c(mNN, mNH, mHN, mHH), 2, 2, byrow=TRUE)
  
#ev <- eigen(M)

#eigenvalues <- ev$values

return(spectralRadius(M))

}

######################################################
# Combined CT (manual on network; digital on both)
######################################################
# assume T_L, T_D are deterministic
# assume D follows a Poisson distribution
######################################################

R_comb = function(mu, beta_N, beta_G, tau, p_D, T_L, T_D, p_M, pi){

#######################################################################################
# Type-1: non-app-user infected without DCT-link, without MCT-link, through network. 
# Type-2: non-app-user infected without DCT-link, with MCT-link, through network.
# Type-3: non-app-user infected without DCT-link, without MCT-link, through homo. 
  
# Type-4: app-user infected without DCT-link, without MCT-link, through network. 
# Type-5: app-user infected without DCT-link, with MCT-link, through network. 
# Type-6: app-user infected with DCT-link,-----, through network. 
# Type-7: app-user infected without DCT-link, without MCT-link, through homo. 
# Type-8: app-user infected with DCT-link, without MCT-link, through homo. 
#######################################################################################
  
D_tilde_minus_1 = mu

m_11 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(1-pi)
m_12 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(p_D*p_M)*(1-pi)  
m_13 = beta_G * tau* (1-pi)
m_14 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(pi)
m_15 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(p_D*p_M)*(pi)
m_16 = 0
m_17 = beta_G * tau* (pi)
m_18 = 0 


if (tau+T_L-T_D <= 0){
  p_r = 1
}
if (tau+T_L-T_D >= tau){
  p_r = 0
}
if (0<tau+T_L-T_D && tau+T_L-T_D < tau){
  p_r = (T_D-T_L)/tau
}

# probability: P(infecting a given friend, traced before recovered)
T_DL = T_D-T_L

if (T_DL < tau && T_DL >= 0){
  p_i_c = (tau-T_DL)/tau - (exp(-beta_N * T_DL)- exp(-beta_N * tau))/(beta_N * tau)
}
if (T_DL < 0 && T_DL > -tau){
  p_i_c = (tau+T_DL)/tau - (1- exp(-beta_N * (tau+T_DL)))/(beta_N * tau)
}
if (T_DL >= tau){
  p_i_c = 0
}
if (T_DL <= -tau){
  p_i_c = 0
}

if (T_DL < tau && T_DL >= 0){
  V = (tau+T_L-T_D)*(tau-T_L+T_D)/(2*tau)
}
if (T_DL < 0 && T_DL > -tau){
  V = (tau-T_L+T_D)^2 /(2*tau)
}
if (T_DL >= tau){
  V = 0
}
if (T_DL <= -tau){
  V = 0
}


m_21 = (D_tilde_minus_1*(1-exp(-beta_N*tau))*p_r*(1-p_M*p_D) + D_tilde_minus_1*p_i_c)*(1-pi)
m_22 = D_tilde_minus_1*(1-exp(-beta_N*tau))*p_D*p_M*p_r*(1-pi)
m_23 = (beta_G * V  + beta_G*tau*p_r)*(1-pi)
m_24 = (D_tilde_minus_1*(1-exp(-beta_N*tau))*p_r*(1-p_M*p_D) + D_tilde_minus_1*p_i_c)*(pi)
m_25 = D_tilde_minus_1*(1-exp(-beta_N*tau))*p_D*p_M*p_r*pi
m_26 = 0 
m_27 = (beta_G * V  + beta_G*tau*p_r)*(pi)
m_28 = 0 


m_31 = mu*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(1-pi)
m_32 = mu*(1-exp(-beta_N*tau))*(p_D*p_M)*(1-pi)
m_33 = beta_G*tau*(1-pi)
m_34 = mu*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(pi)
m_35 = mu*(1-exp(-beta_N*tau))*(p_D*p_M)*(pi)
m_36 = 0
m_37 = beta_G*tau*(pi)
m_38 = 0 


m_41 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(1-pi)
m_42 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(p_D*p_M)*(1-pi)
m_43 = beta_G * tau* (1-pi)
m_44 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-p_D)*(pi)
m_45 = 0
m_46 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(p_D)*(pi)
m_47 = beta_G * tau*(1-p_D)*(pi)
m_48 = beta_G * tau*(p_D)*(pi)

m_51 = m_21
m_52 = m_22
m_53 = m_23
m_54 = (D_tilde_minus_1*(1-exp(-beta_N*tau))*p_r*(1-p_D) + D_tilde_minus_1*p_i_c)*(pi)
m_55 = 0 
m_56 = D_tilde_minus_1*(1-exp(-beta_N*tau))*p_r*(p_D)*pi
m_57 = (beta_G * V  + beta_G*tau*p_r*(1-p_D))*(pi)
m_58 = beta_G*tau*p_r*(p_D)*(pi)


if (T_L < tau){
  p_inf = (tau-T_L)/tau - (1- exp(-beta_N * (tau-T_L)))/(beta_N * tau)
}
if (T_L >= tau){
  p_inf = 0
}

# the average time since being infectious until traced: 

if (T_L < tau){
  T_inf = (tau-T_L)^2 /(2*tau)
}
if (T_L >= tau){
  T_inf = 0
}


m_61 = D_tilde_minus_1*p_inf*(1-pi)
m_62 = 0
m_63 = beta_G * T_inf * (1-pi)
m_64 = D_tilde_minus_1*p_inf*(pi)
m_65 = 0
m_66 = 0 
m_67 = beta_G * T_inf * (pi)
m_68 = 0
  
m_71 = mu*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(1-pi)
m_72 = mu*(1-exp(-beta_N*tau))*(p_D*p_M)*(1-pi)
m_73 = beta_G*tau*(1-pi)
m_74 = mu*(1-exp(-beta_N*tau))*(1-p_D)*pi
m_75 = 0
m_76 = mu*(1-exp(-beta_N*tau))*(p_D)*pi
m_77 = beta_G*tau*(1-p_D)*(pi)
m_78 = beta_G*tau*(p_D)*(pi)

m_81 = mu*p_inf*(1-pi)
m_82 = 0
m_83 = beta_G * T_inf * (1-pi)
m_84 = mu*p_inf*(pi)
m_85 = 0
m_86 = 0 
m_87 =  beta_G * T_inf * (pi)
m_88 = 0



M_combine <- matrix(c(m_11, m_12, m_13, m_14, m_15, m_16, m_17,m_18, 
                      m_21, m_22, m_23, m_24, m_25, m_26,m_27,m_28,
                      m_31, m_32, m_33, m_34, m_35, m_36,m_37,m_38,
                      m_41, m_42, m_43, m_44, m_45, m_46,m_47,m_48,
                      m_51, m_52, m_53, m_54, m_55, m_56,m_57,m_58,
                      m_61, m_62, m_63, m_64, m_65, m_66,m_67,m_68,
                      m_71, m_72, m_73, m_74, m_75, m_76,m_77,m_78,
                      m_81, m_82, m_83, m_84, m_85, m_86,m_87,m_88), 8, 8, byrow=TRUE)
  
#ev <- eigen(M_combine)

#eigenvalues <- ev$values

# Compute the spectral radius of M: the maximum of the absolute values of its eigenvalues
return(spectralRadius(M_combine))
}

######################################################
# Manual with delay on network 
######################################################
# type-1 (N,R): infected through network (N) and reported (R)\\
# type-2 (N,NR): infected through network and not reported (NR)\\
# type-3 (H): infected through homogeneous contacts (H)
######################################################
######################################################

# assume T_L, T_D constant
# assume D follows a Poisson distribution, hence E(D_tilde-1)=mu

R_manual_delay = function(mu, beta_N, beta_G, tau, p_D, T_L, T_D, p_M){
  
D_tilde_minus_1 = mu

############################################  
# elements of the next generation matrix M: 
############################################
  
# For type-1: 

m12 = D_tilde_minus_1*(1-exp(-beta_N*tau))*p_D*p_M
m11 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-p_D*p_M)
m13 = beta_G*tau

# For type-3: 
m32 = mu*(1-exp(-beta_N*tau))*p_D*p_M
m31 = mu*(1-exp(-beta_N*tau))*(1-p_D*p_M)
m33 = beta_G*tau

# For type-2: 
# let 
T_DL = T_D-T_L

# probability: P(recovered before being traced)
if (tau <= T_DL){
  p_r = 1
}
if (T_DL <= 0){
  p_r = 0
}
if (T_DL < tau && 0 < T_DL){
  p_r = (T_D-T_L)/tau
}

# probability: P(infecting a given friend, given that traced before recovered)

if (tau <= T_DL){
  p_i_c = 0
}
if (T_DL <= -tau){
  p_i_c = 0
}
if (T_DL < 0 && T_DL > -tau){
  p_i_c = (tau+T_DL)/tau - (1- exp(-beta_N * (tau+T_DL)))/(beta_N * tau)
}
if (T_DL < tau && T_DL >= 0){
  p_i_c = (tau-T_DL)/tau - (exp(-beta_N * T_DL)- exp(-beta_N * tau))/(beta_N * tau)
}

#print(p_i_c)

# if being traced, the expected infectious period (shortened by CT)

if (tau <= T_DL){
  V = 0
}
if (T_DL <= -tau){
  V = 0
}
if (T_DL < 0 && T_DL > -tau){
  V = (tau-T_L+T_D)^2 /(2*tau)
}
if (T_DL < tau && T_DL >= 0){
  V = (tau+T_L-T_D)*(tau-T_L+T_D)/(2*tau)
}

#print(V)

m21 = D_tilde_minus_1*(1-exp(-beta_N*tau))*p_r*(1-p_M*p_D) + D_tilde_minus_1*p_i_c
m22 = D_tilde_minus_1*(1-exp(-beta_N*tau))*p_r*p_D*p_M
m23 = beta_G * V  + beta_G*tau*p_r
    
M_manual <- matrix(c(m11, m12, m13, m21, m22, m23, m31,m32, m33), 3, 3, byrow=TRUE)
  
#ev <- eigen(M)

#eigenvalues <- ev$values

return(spectralRadius(M_manual))
  
}

######################################################
# Digital CT
######################################################
# pi: fraction of app-users
######################################################
######################################################

# assume T_L, T_D constant
# assume D follows a Poisson distribution

R_digital_network = function(mu, beta_N, beta_G, tau, p_D, T_L, pi){
  
# elements of the next generation matrix M: 

D_tilde_minus_1 = mu
######################################################
# Type-1: non-app-user infected without DCT-link, through network. 
# Type-2: non-app-user infected without DCT-link, through homo. 
# Type-3: app-user infected without DCT-link, through network. 
# Type-4: app-user infected without DCT-link, through homo. 

# Type-5: app-user infected with DCT-link, through network. 
# Type-6: app-user infected with DCT-link, through homo. 
######################################################

m11 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-pi)
m12 = beta_G * tau * (1-pi)
m13 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(pi)
m14 = beta_G * tau * (pi)
m15 = 0
m16 = 0

m21 = mu*(1-exp(-beta_N*tau))*(1-pi)
m22 = beta_G * tau * (1-pi)
m23 = mu*(1-exp(-beta_N*tau))*(pi)
m24 = beta_G * tau * (pi)
m25 = 0
m26 = 0

m31 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-pi)
m32 = beta_G * tau * (1-pi)
m33 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-p_D)*(pi)
m34 = beta_G * tau *(1-p_D) *(pi)
m35 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(p_D)*(pi)
m36 = beta_G * tau *(p_D) *(pi)


m41 = mu*(1-exp(-beta_N*tau))*(1-pi)
m42 = beta_G * tau * (1-pi)
m43 = mu*(1-exp(-beta_N*tau))*(1-p_D)*(pi)
m44 = beta_G * tau *(1-p_D) *(pi)
m45 = mu*(1-exp(-beta_N*tau))*(p_D)*(pi)
m46 = beta_G * tau *(p_D) *(pi)



# probability that a to-be-traced app-user infects a given friend: '
# Note: in this case, the app-user is always traced first before recovered itself

if (T_L < tau){
  p_inf = (tau-T_L)/tau - (1- exp(-beta_N * (tau-T_L)))/(beta_N * tau)
}
if (T_L >= tau){
  p_inf = 0
}

# the average time since being infectious until traced: 

if (T_L < tau){
  T_inf = (tau-T_L)^2 /(2*tau)
}
if (T_L >= tau){
  T_inf = 0
}
 
  
m51 = D_tilde_minus_1*p_inf*(1-pi)
m52 = beta_G * T_inf * (1-pi)
m53 = D_tilde_minus_1*p_inf*(pi)
m54 = beta_G * T_inf *(pi)
m55 = 0
m56 = 0

m61 = mu*p_inf*(1-pi)
m62 = beta_G * T_inf * (1-pi)
m63 = mu*p_inf*(pi)
m64 = beta_G * T_inf *(pi)
m65 = 0
m66 = 0


M_digital <- matrix(c(m11, m12, m13, m14, m15, m16, 
                      m21, m22, m23, m24, m25, m26,
                      m31,m32, m33, m34,m35,m36, 
                      m41,m42,m43,m44,m45,m46,
                      m51,m52,m53,m54,m55,m56,
                      m61,m62,m63,m64,m65,m66), 6, 6, byrow=TRUE)
  
#ev <- eigen(M_digital)

#eigenvalues <- ev$values

# Compute the spectral radius of M: the maximum of the absolute values of its eigenvalues
return(spectralRadius(M_digital))
  
}


######################################################
# Combined CT (manual on network; digital on both)
######################################################
# assume T_L, T_D constant
# assume D follows a Poisson distribution

R_comb = function(mu, beta_N, beta_G, tau, p_D, T_L, T_D, p_M, pi){

######################################################
# Type-1: non-app-user infected without DCT-link, without MCT-link, through network. 
# Type-2: non-app-user infected without DCT-link, with MCT-link, through network.
# Type-3: non-app-user infected without DCT-link, without MCT-link, through homo. 
  
# Type-4: app-user infected without DCT-link, without MCT-link, through network. 
# Type-5: app-user infected without DCT-link, with MCT-link, through network. 
# Type-6: app-user infected with DCT-link,-----, through network. 
# Type-7: app-user infected without DCT-link, without MCT-link, through homo. 
# Type-8: app-user infected with DCT-link, without MCT-link, through homo. 
######################################################
  
D_tilde_minus_1 = mu

m_11 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(1-pi)
m_12 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(p_D*p_M)*(1-pi)  
m_13 = beta_G * tau* (1-pi)
m_14 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(pi)
m_15 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(p_D*p_M)*(pi)
m_16 = 0
m_17 = beta_G * tau* (pi)
m_18 = 0 


if (tau+T_L-T_D <= 0){
  p_r = 1
}
if (tau+T_L-T_D >= tau){
  p_r = 0
}
if (0<tau+T_L-T_D && tau+T_L-T_D < tau){
  p_r = (T_D-T_L)/tau
}

# probability: P(infecting a given friend, traced before recovered)
T_DL = T_D-T_L

if (T_DL < tau && T_DL >= 0){
  p_i_c = (tau-T_DL)/tau - (exp(-beta_N * T_DL)- exp(-beta_N * tau))/(beta_N * tau)
}
if (T_DL < 0 && T_DL > -tau){
  p_i_c = (tau+T_DL)/tau - (1- exp(-beta_N * (tau+T_DL)))/(beta_N * tau)
}
if (T_DL >= tau){
  p_i_c = 0
}
if (T_DL <= -tau){
  p_i_c = 0
}

if (T_DL < tau && T_DL >= 0){
  V = (tau+T_L-T_D)*(tau-T_L+T_D)/(2*tau)
}
if (T_DL < 0 && T_DL > -tau){
  V = (tau-T_L+T_D)^2 /(2*tau)
}
if (T_DL >= tau){
  V = 0
}
if (T_DL <= -tau){
  V = 0
}


m_21 = (D_tilde_minus_1*(1-exp(-beta_N*tau))*p_r*(1-p_M*p_D) + D_tilde_minus_1*p_i_c)*(1-pi)
m_22 = D_tilde_minus_1*(1-exp(-beta_N*tau))*p_D*p_M*p_r*(1-pi)
m_23 = (beta_G * V  + beta_G*tau*p_r)*(1-pi)
m_24 = (D_tilde_minus_1*(1-exp(-beta_N*tau))*p_r*(1-p_M*p_D) + D_tilde_minus_1*p_i_c)*(pi)
m_25 = D_tilde_minus_1*(1-exp(-beta_N*tau))*p_D*p_M*p_r*pi
m_26 = 0 
m_27 = (beta_G * V  + beta_G*tau*p_r)*(pi)
m_28 = 0 


m_31 = mu*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(1-pi)
m_32 = mu*(1-exp(-beta_N*tau))*(p_D*p_M)*(1-pi)
m_33 = beta_G*tau*(1-pi)
m_34 = mu*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(pi)
m_35 = mu*(1-exp(-beta_N*tau))*(p_D*p_M)*(pi)
m_36 = 0
m_37 = beta_G*tau*(pi)
m_38 = 0 


m_41 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(1-pi)
m_42 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(p_D*p_M)*(1-pi)
m_43 = beta_G * tau* (1-pi)
m_44 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(1-p_D)*(pi)
m_45 = 0
m_46 = D_tilde_minus_1*(1-exp(-beta_N*tau))*(p_D)*(pi)
m_47 = beta_G * tau*(1-p_D)*(pi)
m_48 = beta_G * tau*(p_D)*(pi)

m_51 = m_21
m_52 = m_22
m_53 = m_23
m_54 = (D_tilde_minus_1*(1-exp(-beta_N*tau))*p_r*(1-p_D) + D_tilde_minus_1*p_i_c)*(pi)
m_55 = 0 
m_56 = D_tilde_minus_1*(1-exp(-beta_N*tau))*p_r*(p_D)*pi
m_57 = (beta_G * V  + beta_G*tau*p_r*(1-p_D))*(pi)
m_58 = beta_G*tau*p_r*(p_D)*(pi)


if (T_L < tau){
  p_inf = (tau-T_L)/tau - (1- exp(-beta_N * (tau-T_L)))/(beta_N * tau)
}
if (T_L >= tau){
  p_inf = 0
}

# the average time since being infectious until traced: 

if (T_L < tau){
  T_inf = (tau-T_L)^2 /(2*tau)
}
if (T_L >= tau){
  T_inf = 0
}


m_61 = D_tilde_minus_1*p_inf*(1-pi)
m_62 = 0
m_63 = beta_G * T_inf * (1-pi)
m_64 = D_tilde_minus_1*p_inf*(pi)
m_65 = 0
m_66 = 0 
m_67 = beta_G * T_inf * (pi)
m_68 = 0
  
m_71 = mu*(1-exp(-beta_N*tau))*(1-p_D*p_M)*(1-pi)
m_72 = mu*(1-exp(-beta_N*tau))*(p_D*p_M)*(1-pi)
m_73 = beta_G*tau*(1-pi)
m_74 = mu*(1-exp(-beta_N*tau))*(1-p_D)*pi
m_75 = 0
m_76 = mu*(1-exp(-beta_N*tau))*(p_D)*pi
m_77 = beta_G*tau*(1-p_D)*(pi)
m_78 = beta_G*tau*(p_D)*(pi)

m_81 = mu*p_inf*(1-pi)
m_82 = 0
m_83 = beta_G * T_inf * (1-pi)
m_84 = mu*p_inf*(pi)
m_85 = 0
m_86 = 0 
m_87 =  beta_G * T_inf * (pi)
m_88 = 0



M_combine <- matrix(c(m_11, m_12, m_13, m_14, m_15, m_16, m_17,m_18, 
                      m_21, m_22, m_23, m_24, m_25, m_26,m_27,m_28,
                      m_31, m_32, m_33, m_34, m_35, m_36,m_37,m_38,
                      m_41, m_42, m_43, m_44, m_45, m_46,m_47,m_48,
                      m_51, m_52, m_53, m_54, m_55, m_56,m_57,m_58,
                      m_61, m_62, m_63, m_64, m_65, m_66,m_67,m_68,
                      m_71, m_72, m_73, m_74, m_75, m_76,m_77,m_78,
                      m_81, m_82, m_83, m_84, m_85, m_86,m_87,m_88), 8, 8, byrow=TRUE)
  
#ev <- eigen(M_combine)

#eigenvalues <- ev$values

# Compute the spectral radius of M: the maximum of the absolute values of its eigenvalues
return(spectralRadius(M_combine))
  
}

######################################################
# Set epidemic parameters 
######################################################

mu=5
tau=5
T_L=4
################################
# half network transmission:
################################
# we choose beta_G beta_N so that R0 =3

beta_N = log(mu/(mu-0.5*3))/tau
beta_G = (0.5*3)/tau

#other option: we choose beta_G beta_N so that small R0 =1.5

beta_N_s = log(mu/(mu-0.5*1.5))/tau
beta_G_s = (0.5*1.5)/tau

################################
# 75% network transmission:
################################

# we choose beta_G beta_N so that R0 =3

beta_75_N = log(mu/(mu-0.75*3))/tau
beta_75_G = (0.25*3)/tau

#other option: we choose beta_G beta_N so that small R0 =1.5

beta_75_N_s = log(mu/(mu-0.75*1.5))/tau
beta_75_G_s = (0.25*1.5)/tau

# Fix the probability of being diagnosed p_D 
p_D = 0.8

################################################################
# start compute R_DM for each combination of p and pi
################################################################
n_grid = 100

# pi takes values in 
p_app_seq = seq(0,1, length.out=n_grid)
#p_app_seq # increment size with 0.99/99=0.01

# p takes values in 
p_manual_seq = seq(0,1, length.out=n_grid)

# store the R_DM values in a matrix:
R_DM_m=matrix(nrow = n_grid,ncol = n_grid)
R_DM_s_m=matrix(nrow = n_grid,ncol = n_grid)
R_DM_75m=matrix(nrow = n_grid,ncol = n_grid)
R_DM_s_75m=matrix(nrow = n_grid,ncol = n_grid)

for (k in (1:n_grid)) {
  
  p_A = p_app_seq[k]
  
  for (j in (1:n_grid)) {
    
  p = p_manual_seq[j] 
  
  R_DM_m[j,k] = R_comb(mu=mu, beta_N=beta_N, beta_G= beta_G, 
             tau=tau, p_D=p_D, T_L= T_L, T_D=3,
             p=p,pi=p_A)
  
  R_DM_75m[j,k] = R_comb(mu=mu, beta_N=beta_75_N, beta_G= beta_75_G, 
             tau=tau, p_D=p_D, T_L= T_L, T_D=3,
             p=p,pi=p_A)
  
  R_DM_s_m[j,k] = R_comb(mu=mu, beta_N=beta_N_s, beta_G= beta_G_s, 
             tau=tau, p_D=p_D, T_L= T_L, T_D=3,
             p=p,pi=p_A)
  
  R_DM_s_75m[j,k] = R_comb(mu=mu, beta_N=beta_75_N_s, beta_G= beta_75_G_s, 
             tau=tau, p_D=p_D, T_L= T_L, T_D=3,
             p=p,pi=p_A)
  
  }
}

################################
# store in a data frame 
################################

data_comb <- expand.grid(manual_prob=p_manual_seq,app_prob=p_app_seq)
data_comb$R_MD = c(R_DM_m)
data_comb$R_MD_s = c(R_DM_s_m)
data_comb$R_MD_75 = c(R_DM_75m)
data_comb$R_MD_s_75 = c( R_DM_s_75m)


################################
# compute R0(1-r_D)(1-r_M)
################################

R_naive_m=matrix(nrow = n_grid,ncol = n_grid)
R_naive_m_75=matrix(nrow = n_grid,ncol = n_grid)
R_naive_ms=matrix(nrow = n_grid,ncol = n_grid)
R_naive_ms_75=matrix(nrow = n_grid,ncol = n_grid)


for (k in (1:(n_grid))) {
  
  p_A = p_app_seq[k]
  
  for (j in (1:(n_grid))) {
  p = p_manual_seq[j]  
  
  # when R0 = 3, 50% transmission on network
  
  reduction_mct = R_manual_delay(mu=mu, beta_N= beta_N, beta_G= beta_G, tau=tau, p=p, p_D=0.8, T_L= T_L, T_D=3)/3
  
  reduction_dct =R_digital_network(mu=mu, beta_N= beta_N, beta_G= beta_G, tau=tau, pi=p_A, p_D=0.8, T_L= T_L)/3
  
  # when R0 = 1.5, 50% transmission on network
  
  reduction_mct_s = R_manual_delay(mu=mu, beta_N= beta_N_s, beta_G= beta_G_s, tau=tau, p=p, p_D=0.8, T_L= T_L, T_D=3)/1.5
  
  reduction_dct_s = R_digital_network(mu=mu, beta_N= beta_N_s, beta_G= beta_G_s, tau=tau, pi=p_A, p_D=0.8, T_L= T_L)/1.5
  
  # when R0 = 3, 75% transmission on network
  
  reduction_mct_75 = R_manual_delay(mu=mu, beta_N= beta_75_N, beta_G= beta_75_G, tau=tau, p=p, p_D=0.8, T_L= T_L, T_D=3)/3
  
  reduction_dct_75 =R_digital_network(mu=mu, beta_N= beta_75_N, beta_G= beta_75_G, tau=tau, pi=p_A, p_D=0.8, T_L= T_L)/3
  
  # when R0 = 1.5, 75% transmission on network
  
  reduction_mct_s_75 = R_manual_delay(mu=mu, beta_N= beta_75_N_s, beta_G= beta_75_G_s, tau=tau, p=p, p_D=0.8, T_L= T_L, T_D=3)/1.5
  
  reduction_dct_s_75 = R_digital_network(mu=mu, beta_N= beta_75_N_s, beta_G= beta_75_G_s, tau=tau, pi=p_A, p_D=0.8, T_L= T_L)/1.5
  
  
  
  R_naive_m[j,k] = 3*reduction_mct*reduction_dct
  
  R_naive_ms[j,k] = 1.5*reduction_mct_s*reduction_dct_s
  
  R_naive_m_75[j,k] = 3*reduction_mct_75*reduction_dct_75
  
  R_naive_ms_75[j,k] = 1.5*reduction_mct_s_75*reduction_dct_s_75
  
  }
  
}

################################
# store in the data frame 
################################

data_comb$R_naive <-c(R_naive_m)
data_comb$R_naive_s <-c(R_naive_ms)
data_comb$R_naive_75 <-c(R_naive_m_75)
data_comb$R_naive_s_75 <-c(R_naive_ms_75)

#####################################################
# compare critical lines: R_DM=1 with R0(1-r_D)(1-r_M)=1
#####################################################

plot_R_MD_naive = ggplot(data_comb)+ 
labs(title="", y="app-using fraction pi_A", x="manual reporting probability p_M")+
# set contour breaks at desired level R_DM =1
geom_contour(aes(y=app_prob,x=manual_prob,z = R_MD), 
               breaks = 1,size = 0.5,col= 'hotpink')+
# set contour breaks at desired level R_toy =1
geom_contour(aes(y=app_prob,x=manual_prob,z = R_naive), 
               breaks = 1, col = 'hotpink',size = 0.5, linetype ="dashed")+ 
# set contour breaks at desired level R_DM =1
geom_contour(aes(y=app_prob,x=manual_prob,z = R_MD_75), 
               breaks = 1,size = 0.5,col= 'dodgerblue')+
# set contour breaks at desired level R_toy =1
geom_contour(aes(y=app_prob,x=manual_prob,z = R_naive_75), 
               breaks = 1, col = 'dodgerblue',size = 0.5, linetype ="dashed")+ 
# set contour breaks at desired level R_DM =1
geom_contour(aes(y=app_prob,x=manual_prob,z = R_MD_s), 
               breaks = 1,size = 0.5,col= 'darkorange')+
# set contour breaks at desired level R_toy =1
geom_contour(aes(y=app_prob,x=manual_prob,z = R_naive_s), 
               breaks = 1, col = 'darkorange',size = 0.5, linetype ="dashed")+ 
# set contour breaks at desired level R_DM =1
geom_contour(aes(y=app_prob,x=manual_prob,z = R_MD_s_75), 
               breaks = 1,size = 0.5,col= 'firebrick1')+
# set contour breaks at desired level R_toy =1
geom_contour(aes(y=app_prob,x=manual_prob,z = R_naive_s_75), 
               breaks = 1, col = 'firebrick1',size = 0.5, linetype ="dashed")+ 
theme(plot.title = element_text(hjust = 0.5,face="bold"))+
ylim(0, 1)+
theme(text = element_text(size=14))

#####################################################
# produce a legend
#####################################################

# Create a DataFrame 
df_sample <- data.frame(
  Xdata = rnorm(4),                        
  Ydata = rnorm(4),
  LegendData = c("R0=3; alpha=50%", 
                 "R0=3; alpha=75%",  
                 "R0=1.5; alpha=50%", 
                 "R0=1.5; alpha=75%"))
  
# Create a Scatter Plot and assign it 
# to gplot data object
gplot <- ggplot(df_sample, aes(Xdata, Ydata, color = LegendData)) +   
  geom_line()+
  scale_color_manual(values=c('hotpink','dodgerblue', 'darkorange','firebrick1'), 
                       name="",
                       labels=c("R0=3; alpha=50%", 
                 "R0=3; alpha=75%",  
                 "R0=1.5; alpha=50%", 
                 "R0=1.5; alpha=75%"))+
  theme(text = element_text(size=14))


legend = get_legend(gplot)

#####################################################
# add this legend to the plot
#####################################################

grid.arrange(plot_R_MD_naive,legend,ncol=2, widths=c(3, 0.8))