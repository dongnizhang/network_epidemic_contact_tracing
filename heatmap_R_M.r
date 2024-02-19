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

# Set epidemic parameters 

mu=5
tau=5
T_L=4

# half network transmission:
# we choose beta_G beta_N so that R0 =3

beta_N = log(mu/(mu-0.5*3))/tau
beta_G = (0.5*3)/tau

#other option: we choose beta_G beta_N so that small R0 =1.5

beta_N_s = log(mu/(mu-0.5*1.5))/tau
beta_G_s = (0.5*1.5)/tau

# we choose beta_G beta_N so that R0 =3

beta_75_N = log(mu/(mu-0.75*3))/tau
beta_75_G = (0.25*3)/tau

#other option: we choose beta_G beta_N so that small R0 =1.5

beta_75_N_s = log(mu/(mu-0.75*1.5))/tau
beta_75_G_s = (0.25*1.5)/tau



# 50% transmission on network: Heatmap of R_M (fix T_D, vary p_M and p_D)


n_grid=100

# Fix delay to be 3 days
T_D = 3
#T_D = 1

# p_M takes values in 
p_seq = seq(0,1, length.out=n_grid)

# p_D takes values in 
p_D_seq = seq(0.01,1, length.out=n_grid)

# store the R_M values in a matrix:
R_M_m=matrix(nrow = n_grid,ncol = n_grid)

# store the R_M values in a matrix (for small R0)
R_M_s_m=matrix(nrow = n_grid,ncol = n_grid)


for (k in (1:n_grid)) {
  
  p1 = p_seq[k]
  
  for (j in (1:n_grid)) {
    
    p2 = p_D_seq[j] 
    
    R = R_manual_delay(mu=mu, beta_N= beta_N, beta_G= beta_G, 
                       tau=tau, T_L= T_L, T_D=T_D,
                       p_M=p1,p_D=p2)
    
    
    R_M_m[j,k] = R
    
    Rs = R_manual_delay(mu=mu, beta_N= beta_N_s, beta_G= beta_G_s, 
                        tau=tau, T_L= T_L, T_D=T_D,
                        p_M=p1,p_D=p2)
    R_M_s_m[j,k] = Rs
    
    
  }
}

data_fix_delay <- expand.grid(diag_prob=p_D_seq,prob=p_seq)
data_fix_delay$R_M = c(R_M_m)
data_fix_delay$R_M_s = c(R_M_s_m)


# 75% transmission on network: Heatmap of R_M (fix T_D, vary p_M and p_D)


n_grid=100

# Fix delay to be 3 days
T_D = 3
#T_D=1

# p_M takes values in 
p_seq = seq(0,1, length.out=n_grid)

# p_D takes values in 
p_D_seq = seq(0.01,1, length.out=n_grid)

# store the R_M values in a matrix:
R_M_75_m=matrix(nrow = n_grid,ncol = n_grid)

# store the R_M values in a matrix (for small R0)
R_M_75_s_m=matrix(nrow = n_grid,ncol = n_grid)


for (k in (1:n_grid)) {
  
  p1 = p_seq[k]
  
  for (j in (1:n_grid)) {
    
    p2 = p_D_seq[j] 
    
    R = R_manual_delay(mu=mu, beta_N= beta_75_N, beta_G= beta_75_G, 
                       tau=tau, T_L= T_L, T_D=T_D,
                       p_M=p1,p_D=p2)
    
    
    R_M_75_m[j,k] = R
    
    Rs = R_manual_delay(mu=mu, beta_N= beta_75_N_s, beta_G= beta_75_G_s, 
                        tau=tau, T_L= T_L, T_D=T_D,
                        p_M=p1,p_D=p2)
    R_M_75_s_m[j,k] = Rs
    
    
  }
}

data_fix_delay$R_M_75 = c(R_M_75_m)
data_fix_delay$R_M_75_s = c(R_M_75_s_m)



plot_p_D_p_M_R_M = ggplot(data_fix_delay, aes(y=diag_prob,x=prob,fill=R_M))+ 
  geom_tile()+
  labs(title="Heatmap of R_M when R0=3", subtitle = "(50% network transmission)", y="probability of diagnosis p_D", x="manual reporting probability p_M")+
  # set contour breaks at desired level
  geom_contour(aes(z = R_M), breaks = 2.25, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M), breaks = 2.75, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M), 
               breaks = 2, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M), 
               breaks = 2.5, col = 'black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_M), stroke = 0.1,breaks=2.25,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M), stroke = 0.1,breaks=2.75,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M), stroke = 0.1,breaks=2.5,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M), stroke = 0.1,breaks=2,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5,face = "bold"))+
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,R_0(mu,beta_N = beta_N, beta_G = beta_G, tau=tau)),"R_M")+
  ylim(0.01, 1)

plot_p_D_p_M_R_M_s = ggplot(data_fix_delay, aes(y=diag_prob,x=prob,fill=R_M_s))+ 
  geom_tile()+
  labs(title="Heatmap of R_M when R0=1.5", subtitle = "(50% network transmission)", y="probability of diagnosis p_D", x="manual reporting probability p_M")+
  # set contour breaks at desired level
  geom_contour(aes(z = R_M_s), breaks = 1, col = 'black',size = 0.5)+
  geom_contour(aes(z = R_M_s), breaks = 1.4, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_s), breaks = 1.3, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_s), breaks = 1.2, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_s), breaks = 1.1, col = 'black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_M_s), stroke = 0.1,breaks=1,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_s), stroke = 0.1,breaks=1.1,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_s), stroke = 0.1,breaks=1.2,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_s), stroke = 0.1,breaks=1.3,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_s), stroke = 0.1,breaks=1.4,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5,face = "bold"))+
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,R_0(mu,beta_N = beta_N, beta_G = beta_G, tau=tau)),"R_M")+
  ylim(0.01, 1)

plot_p_D_p_M_R_M_75 = ggplot(data_fix_delay, aes(y=diag_prob,x=prob,fill=R_M_75))+ 
  geom_tile()+
  labs(title="Heatmap of R_M when R0=3", subtitle="(75% network transmission)", y="probability of diagnosis p_D", x="manual reporting probability p_M")+
  # set contour breaks at desired level
  geom_contour(aes(z = R_M_75), breaks = 2.75, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75), breaks = 2.25, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75), breaks = 1.75, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75), 
               breaks = 2, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75), 
               breaks = 2.5, col = 'black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=1.75,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=2.25,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=2.75,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=2.5,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=2,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5,face = "bold"))+
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,R_0(mu,beta_N = beta_N, beta_G = beta_G, tau=tau)),"R_M")+
  ylim(0.01, 1)

plot_p_D_p_M_R_M_75s = ggplot(data_fix_delay, aes(y=diag_prob,x=prob,fill=R_M_75_s))+ 
  geom_tile()+
  labs(title="Heatmap of R_M when R0=1.5",subtitle="(75% network transmission)", y="probability of diagnosis p_D", x="manual reporting probability p_M", )+
  # set contour breaks at desired level
  geom_contour(aes(z = R_M_75_s), breaks = 1, col = 'black',size = 0.5)+
  geom_contour(aes(z = R_M_75_s), breaks = 1.1, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75_s), breaks = 1.2, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75_s), breaks = 1.3, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75_s), breaks = 1.4, col = 'black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_M_75_s), stroke = 0.1,breaks=1,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_75_s), stroke = 0.1,breaks=1.1,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_75_s), stroke = 0.1,breaks=1.2,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_75_s), stroke = 0.1,breaks=1.3,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_75_s), stroke = 0.1,breaks=1.4,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5,face = "bold"))+
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,R_0(mu,beta_N = beta_N, beta_G = beta_G, tau=tau)),"R_M")+
  ylim(0.01, 1)




# Heatmap of R_M (fix p_D, vary T_D and p_M)

n_grid=100

# fix p_D = 0.8
p_D = 0.8

# p_M takes values in 
p_seq = seq(0,1, length.out=n_grid)

# T_D takes values in 
T_D_seq = seq(0,6, length.out=n_grid)

# store the R_M values in a matrix:
R_M_m=matrix(nrow = n_grid,ncol = n_grid)

# store the R_M values in a matrix (for small R0)
R_M_s_m=matrix(nrow = n_grid,ncol = n_grid)

# store the R_M values in a matrix:
R_M_75_m=matrix(nrow = n_grid,ncol = n_grid)

# store the R_M values in a matrix (for small R0)
R_M_75_s_m=matrix(nrow = n_grid,ncol = n_grid)


for (k in (1:n_grid)) {
  
  p = p_seq[k]
  
  for (j in (1:n_grid)) {
    
    t = T_D_seq[j] 
    
    R = R_manual_delay(mu=mu, beta_N= beta_N, beta_G= beta_G, 
                       tau=tau, T_L= T_L, p_D=p_D,
                       p_M=p,T_D=t)
    
    R2 = R_manual_delay(mu=mu, beta_N= beta_75_N, beta_G= beta_75_G, 
                        tau=tau, T_L= T_L, p_D=p_D,
                        p_M=p,T_D=t)
    
    R_M_m[j,k] = R
    
    R_M_75_m[j,k] = R2
    
    
    Rs = R_manual_delay(mu=mu, beta_N= beta_N_s, beta_G= beta_G_s, 
                        tau=tau, T_L= T_L, p_D=p_D,
                        p_M=p,T_D=t)
    Rs2 = R_manual_delay(mu=mu, beta_N= beta_75_N_s, beta_G= beta_75_G_s, 
                         tau=tau, T_L= T_L, p_D=p_D,
                         p_M=p,T_D=t)
    
    R_M_s_m[j,k] = Rs
    
    R_M_75_s_m[j,k] = Rs2
    
  }
}

data_fix_p_D <- expand.grid(delay=T_D_seq,prob=p_seq)
data_fix_p_D$R_M = c(R_M_m)
data_fix_p_D$R_M_s = c(R_M_s_m)
data_fix_p_D$R_M_75 = c(R_M_75_m)
data_fix_p_D$R_M_s_75 = c(R_M_75_s_m)



plot_T_D_p_M_R_M = ggplot(data_fix_p_D, aes(y=delay,x=prob,fill=R_M))+ 
  geom_tile()+
  labs(title="Heatmap of R_M when R0=3", 
       subtitle = "(50% network transmission)", 
       y="tracing delay", x="manual reporting probability p_M")+
  # set contour breaks at desired level
  geom_contour(aes(z = R_M), 
               breaks = 2.75, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M), 
               breaks = 2.5, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M), 
               breaks = 2.25, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M), 
               breaks = 2, col = 'black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_M), stroke = 0.1,breaks=2.75,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M), stroke = 0.1,breaks=2.5,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M), stroke = 0.1,breaks=2.25,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M), stroke = 0.1,breaks=2,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5,face = "bold"))+
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,R_0(mu,beta_N = beta_N, beta_G = beta_G, tau=tau)),"R_M")+
  ylim(0,6)

plot_T_D_p_M_R_M_s = ggplot(data_fix_p_D, aes(y=delay,x=prob,fill=R_M_s))+ 
  geom_tile()+
  labs(title="Heatmap of R_M when R0=1.5", 
       subtitle = "(50% network transmission)",  
       y="tracing delay", x="manual reporting probability p_M")+
  # set contour breaks at desired level
  geom_contour(aes(z = R_M_s), breaks = 1, col = 'black',size = 0.5)+
  geom_contour(aes(z = R_M_s), breaks = 1.1, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_s), breaks = 1.2, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_s), breaks = 1.3, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_s), breaks = 1.4, col = 'black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_M_s), stroke = 0.1,breaks=1.1,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_s), stroke = 0.1,breaks=1.2,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_s), stroke = 0.1,breaks=1.3,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_s), stroke = 0.1,breaks=1.4,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_s), stroke = 0.1,breaks=1,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5,face = "bold"))+
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,R_0(mu,beta_N = beta_N, beta_G = beta_G, tau=tau)),"R_M")+
  ylim(0,6)

plot_T_D_p_M_R_M_75 = ggplot(data_fix_p_D, aes(y=delay,x=prob,fill=R_M_75))+ 
  geom_tile()+
  labs(title="Heatmap of R_M when R0=3", subtitle = "(75% network transmission)", y="tracing delay", x="manual reporting probability p_M")+
  # set contour breaks at desired level
  geom_contour(aes(z = R_M_75), breaks = 1.5, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75), breaks = 1.75, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75), breaks = 2, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75), breaks = 2.25, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75), breaks = 2.5, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_75), breaks = 2.75, col = 'black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=1.5,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=1.75,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=2,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=2.25,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=2.5,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_75), stroke = 0.1,breaks=2.75,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5,face = "bold"),text = element_text(size=10))+
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,R_0(mu,beta_N = beta_N, beta_G = beta_G, tau=tau)),"R_M")+
  ylim(0,6)



plot_T_D_p_M_R_M_75s = ggplot(data_fix_p_D, aes(y=delay,x=prob,fill=R_M_s_75))+ 
  geom_tile()+
  labs(title="Heatmap of R_M when R0=1.5",subtitle="(75% network transmission)", y="tracing delay", x="manual reporting probability p_M")+
  # set contour breaks at desired level
  geom_contour(aes(z = R_M_s_75), breaks = 1, col = 'black',size = 0.5)+
  geom_contour(aes(z = R_M_s_75), breaks = 1.1, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_s_75), breaks = 1.2, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_s_75), breaks = 1.3, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_M_s_75), breaks = 1.4, col = 'black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_M_s_75), stroke = 0.1,breaks=1.1,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_s_75), stroke = 0.1,breaks=1.2,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_s_75), stroke = 0.1,breaks=1.3,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  geom_text_contour(aes(z = R_M_s_75), stroke = 0.1,breaks=1.4,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_M_s_75), stroke = 0.1,breaks=1,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold"),plot.subtitle = element_text(hjust = 0.5,face = "bold"))+
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,R_0(mu,beta_N = beta_N, beta_G = beta_G, tau=tau)),"R_M")+
  ylim(0,6)



# save the four heatmaps (R0=3,50/75% transmission)

# 2. Save the legend
#+++++++++++++++++++++++
legend = get_legend(plot_T_D_p_M_R_M)
# 3. Remove the legend 
#+++++++++++++++++++++++
plot_T_D_p_M_R_M <- plot_T_D_p_M_R_M + theme(legend.position="none")
plot_T_D_p_M_R_M_75 <- plot_T_D_p_M_R_M_75 + theme(legend.position="none")
plot_p_D_p_M_R_M <- plot_p_D_p_M_R_M + theme(legend.position="none")
plot_p_D_p_M_R_M_75 <- plot_p_D_p_M_R_M_75 + theme(legend.position="none")


# 4. Arrange ggplot2 graphs with a specific width
grid.arrange(plot_p_D_p_M_R_M,plot_p_D_p_M_R_M_75, legend,plot_T_D_p_M_R_M, plot_T_D_p_M_R_M_75,ncol=3, widths=c(2.3, 2.3,0.8))



# save the four heatmaps (R0=1.5,50/75% transmission)

#  Remove the legend 
#+++++++++++++++++++++++
plot_T_D_p_M_R_M_s <- plot_T_D_p_M_R_M_s + theme(legend.position="none")
plot_T_D_p_M_R_M_75s <- plot_T_D_p_M_R_M_75s + theme(legend.position="none")

plot_p_D_p_M_R_M_s <- plot_p_D_p_M_R_M_s + theme(legend.position="none")
plot_p_D_p_M_R_M_75s <- plot_p_D_p_M_R_M_75s + theme(legend.position="none")


# 4. Arrange ggplot2 graphs with a specific width
grid.arrange(plot_p_D_p_M_R_M_s,plot_p_D_p_M_R_M_75s, legend,plot_T_D_p_M_R_M_s, plot_T_D_p_M_R_M_75s,ncol=3, widths=c(2.3, 2.3,0.8))



