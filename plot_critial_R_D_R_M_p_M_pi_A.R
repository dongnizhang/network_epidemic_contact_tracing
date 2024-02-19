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

# Basic reproduction number: without CT

R_0 = function(mu, beta_N, beta_G, tau){
  
  ## assume D follows a poisson distribution.
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

# Manual Reproduction Number

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

# Digital reproduction number

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

# Set epidemic parameters 

mu=5
tau=5
T_L=4

# half network transmission:
# we choose beta_G beta_N so that R0 =3

beta_N = log(mu/(mu-0.5*3))/tau
beta_G = (0.5*3)/tau


# plot R_D and R_M against pi, pi^2 and p: 

pi_vector = seq(0.01,1, length.out=100)
p_vector  = seq(0.01,1, length.out=100)
r_D_vec = c()
r_M_vec = c()
r_D_sq_vec = c()

for (k in 1:length(pi_vector)){
  
  x = R_digital_network(mu=5, beta_N=beta_N, beta_G=beta_G, tau=5, pi=pi_vector[k], p_D=0.8, T_L=4)
  
  
  y = R_manual_delay(mu=5, beta_N=beta_N, beta_G=beta_G, tau=5, p=p_vector[k], p_D=0.8, T_L=4, T_D=3)
  
  z = R_digital_network(mu=5, beta_N=beta_N, beta_G=beta_G, tau=5, pi=sqrt(pi_vector[k]), p_D=0.8, T_L=4)
  
  r_D_vec = c(r_D_vec,(3-x)/3)
  r_D_sq_vec = c(r_D_sq_vec,(3-z)/3)
  r_M_vec = c(r_M_vec,(3-y)/3)
}

data_compare = data.frame(
  prob = c(p_vector,pi_vector,pi_vector),
  r = c(r_M_vec,r_D_vec,r_D_sq_vec),
  r_type = c(rep("r_M against p_M", 100), rep("r_D against pi_A", 100), rep("r_D against (pi_A)^2", 100)))

ggplot(data_compare, aes(x = prob, y = r, group = r_type, col = r_type,linetype = r_type))+ 
  geom_line()+
  ylab("r_M, r_D")+
  xlab("p_M, pi_A, (pi_A)^2")+
  scale_color_manual("", values = c("r_M against p_M"="#ff9200", 
                                    "r_D against pi_A"="#0094ff", 
                                    "r_D against (pi_A)^2"="#0094ff"))+
  scale_linetype_manual("", values = c("r_M against p_M"="solid", 
                                       "r_D against pi_A"="solid", 
                                       "r_D against (pi_A)^2"="dashed"))+
  xlim(0.01,1)+
  theme(text = element_text(size=14))