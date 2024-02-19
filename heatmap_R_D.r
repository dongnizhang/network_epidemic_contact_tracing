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

#other option: we choose beta_G beta_N so that small R0 =1.5

beta_N_s = log(mu/(mu-0.5*1.5))/tau
beta_G_s = (0.5*1.5)/tau

# we choose beta_G beta_N so that R0 =3

beta_75_N = log(mu/(mu-0.75*3))/tau
beta_75_G = (0.25*3)/tau

#other option: we choose beta_G beta_N so that small R0 =1.5

beta_75_N_s = log(mu/(mu-0.75*1.5))/tau
beta_75_G_s = (0.25*1.5)/tau




# Heat map of R_D for each combination of pi and p_D

n_grid = 100

# pi takes values in 
pi_seq = seq(0,1, length.out=n_grid)

# p_D takes values in 
p_D_seq = seq(0.01,1, length.out=n_grid)

# store the R_D values in a matrix:
R_D_m=matrix(nrow = n_grid,ncol = n_grid)
# store the R_D values in a matrix (small R0)
R_D_s_m=matrix(nrow = n_grid,ncol = n_grid)

for (k in (1:n_grid)) {
  
  p1 = pi_seq[k]
  
  for (j in (1:n_grid)) {
    
    p2 = p_D_seq[j] 
    
    R = R_digital_network(mu=mu, beta_N= beta_N, beta_G= beta_G, 
                          tau=tau, p_D=p2, T_L= T_L,
                          pi=p1)
    
    Rs = R_digital_network(mu=mu, beta_N= beta_N_s, beta_G= beta_G_s, 
                           tau=tau, p_D=p2, T_L= T_L,
                           pi=p1)
    
    R_D_m[j,k] = R
    R_D_s_m[j,k] = Rs
    
  }
}

# store all in a data frame 
data_dct <- expand.grid(diag_prob=p_D_seq,prob=pi_seq)
data_dct$R_D = c(R_D_m)
data_dct$R_D_s = c(R_D_s_m)


# Heatmap for R_D (vary with p_D and pi)


plot_R_D = ggplot(data_dct, aes(y=diag_prob,x=prob,fill=R_D))+ 
  geom_tile()+
  labs(title="Heatmap of R_D when R0 = 3", y="probability of diagnosis p_D", x="app-using fraction pi_A")+
  # set contour breaks at desired level
  geom_contour(aes(z = R_D), breaks = 1, col ='black',size = 0.5)+
  geom_contour(aes(z = R_D), 
               breaks = 1.5, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_D), 
               breaks = 2, col = 'black',size = 0.5,linetype ="dashed")+
  geom_contour(aes(z = R_D), 
               breaks = 2.5, col = 'black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_D), stroke = 0.1,breaks=2.5,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_D), stroke = 0.1,breaks=2,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_D), stroke = 0.1,breaks=1.5,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_text_contour(aes(z = R_D), stroke = 0.1,breaks=1,
                    label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,R_0(mu,beta_N = beta_N, beta_G = beta_G, tau=tau)),"R_D")+
  ylim(0.01, 1)


# heatmap for Reproduction Number for digital tracing (smaller R0)
plot_R_D_s = ggplot(data_dct, aes(y=diag_prob,x=prob,fill=R_D_s))+ 
  geom_tile()+
  labs(title="Heatmap of R_D when R0 = 1.5", y="probability of diagnosis p_D", x="app-using fraction pi_A")+
  # set contour breaks at desired level
  geom_contour(aes(z = R_D_s), breaks = 1, col ='black',size = 0.5)+
  geom_text_contour(aes(z = R_D_s), stroke = 0.1,breaks=1,label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_contour(aes(z = R_D_s), breaks = 1.4, col ='black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_D_s), breaks = 1.4, stroke = 0.1,label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_contour(aes(z = R_D_s), breaks = 1.3, col ='black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_D_s), breaks = 1.3, stroke = 0.1,label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_contour(aes(z = R_D_s), breaks = 1.2, col ='black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_D_s), breaks = 1.2, stroke = 0.1,label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  geom_contour(aes(z = R_D_s), breaks = 1.1, col ='black',size = 0.5,linetype ="dashed")+
  geom_text_contour(aes(z = R_D_s), breaks = 1.1, stroke = 0.1,label.placer=label_placer_fraction(frac = 0.5, rot_adjuster = isoband::angle_halfcircle_bottom()))+
  theme(text = element_text(size=10))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_distiller(palette = "RdYlBu",limits = c(0,R_0(mu,beta_N = beta_N, beta_G = beta_G, tau=tau)),"R_D")+
  ylim(0.01, 1)

# 2. Save the legend
#+++++++++++++++++++++++
legend = get_legend(plot_R_D)
# 3. Remove the legend 
#+++++++++++++++++++++++
plot_R_D <- plot_R_D + theme(legend.position="none")
plot_R_D_s <- plot_R_D_s + theme(legend.position="none")

# 4. Arrange ggplot2 graphs with a specific width
grid.arrange(plot_R_D, plot_R_D_s,legend,ncol=3, widths=c(2.3, 2.3,0.8))





