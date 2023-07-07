setwd("~/OneDrive - The University of Nottingham/Scripts/PairwiseDSA")
rm(list=ls())

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(ggplot2); 
library(deSolve)

source("tudColours.R")

############################################ ODE for SIR 
pairwise.ODE <- function(Tmax, params, ic, dt=0.01){
  R0 <- params["R0"]
  n <- params["n"]
  gamma <- params["gamma"]
  rho <- params["rho"]
  xi = 1.0*(n-1)/n;
  tau = R0*gamma/(n*xi - R0);
  b = n*tau; 
  a = tau + gamma;
  
  p=length(ic)
  n=Tmax/dt
  xmat = matrix(0,ncol=(p+1), nrow=(n+1))
  x=ic
  xmat[1,]=c(0,x)
  for(i in 2:(n+1)){
    x = x+c(- a*x[1]*(1-x[1]**(xi-1))/(1-xi) - b*(1- (x[1]**xi))*(x[1]**xi) - b*rho*(x[1]**xi))*dt
    xmat[i,2:(p+1)]=x
    xmat[i,1]=xmat[(i-1),1]+dt
  }
  return(xmat)
}

library(deSolve)

S_dot <- function(t, x, params){
  S <- x[1]
  R0 <- params["R0"]
  n <- params["n"]
  gamma <- params["gamma"]
  rho <- params["rho"]
  xi = 1.0*(n-1)/n;
  tau = R0*gamma/(n*xi - R0);
  b = n*tau; 
  a = tau + gamma;
  dSdt <- - a*x[1]*(1-x[1]**(xi-1))/(1-xi) - b*(1- (x[1]**xi))*(x[1]**xi) - b*rho*(x[1]**xi)
  dxdt <- c(dSdt)
  list(dxdt)
}


#### Other inputs;
population_size <- 10000
# I0 <- 100
I0 <- 1

## Stan configurations
N <- 5000 #sample size for infection data
M <- 5000 #sample size for recovery data 
nIter <- 5000 #number of iterations of Hamiltonian Monte Carlo 
nChains <- 2 #number of parallel chains

sellke.data <-read.csv("data/epidata_netw_1_run_1.csv", header =  FALSE)
colnames(sellke.data) <- c("Time", "S", "I")
# colnames(sellke.data) <- c("Time", "S", "I", "Vertex")


L <- length(sellke.data[,1])
for (i in 2:L){
  if (sellke.data$S[i]<sellke.data$S[i-1]){
    sellke.data$EventType[i] = 'Infection'
  } else {
    sellke.data$EventType[i] = 'Recovery'
  }
}

all_infection_times <- sellke.data[which(sellke.data$S < population_size - I0+1 & sellke.data$S > I0-1 &  sellke.data$EventType=="Infection"), 1, drop = FALSE]$Time
# all_infection_times <- sellke.data[which(sellke.data$EventType=="Infection"), 1, drop = FALSE]$Time

all_infection_times <- all_infection_times[is.na(all_infection_times) == FALSE]

initial_time <- min(all_infection_times)
# initial_time <- 0
all_infection_times <- all_infection_times - initial_time
all_infection_times <- all_infection_times[all_infection_times >0]

all_recovery_times <- sellke.data[which(sellke.data$Time> initial_time &  sellke.data$EventType=="Recovery"), 1, drop = FALSE]$Time
all_recovery_times <- all_recovery_times[is.na(all_recovery_times) == FALSE]
all_recovery_times <- all_recovery_times - initial_time

### Take a random sample of infection and recovery times
if (length(all_infection_times) < N){
  infection_times_sample <- sort(unique(all_infection_times), decreasing = FALSE)
  N <- length(infection_times_sample)
} else {
  infection_times_sample <- sort(sample(unique(all_infection_times), size = N, replace = FALSE), decreasing = FALSE)
  #infection_times_sample <- sort(unique(infection_times_sample), decreasing = FALSE)
  N <- length(infection_times_sample)
}

if (length(all_recovery_times) < M){
  recovery_times_sample <- sort(unique(all_recovery_times), decreasing = FALSE)
  M <- length(recovery_times_sample)
} else {
  recovery_times_sample <- sort(sample(unique(all_recovery_times), size = M, replace = FALSE), decreasing = FALSE)
  #infection_times_sample <- sort(unique(infection_times_sample), decreasing = FALSE)
  N <- length(recovery_times_sample)
}


hist(infection_times_sample)

hist(recovery_times_sample)

#### Create Stan data

stan_data <- list(N = N, 
                  infection_times = infection_times_sample,
                  M = M,
                  recovery_times = recovery_times_sample)

file <- file.path(getwd(),"pairwise_uninformative.stan")
mod <- cmdstan_model(file)

print(mod$print())

fit <- mod$sample(
  data = stan_data, 
  seed = 123, 
  chains = nChains, 
  parallel_chains = nChains,
  refresh = 500 # print update every 500 iters
)

fit$summary()

print(fit$summary())

posterior_samples <- fit$draws(format = "df")



# posterior_samples <- read.csv("plots/synthetic/synthetic_DSA_posterior_samples.csv")


R0_plot<-ggplot(posterior_samples, aes(x=R0)) + 
  geom_histogram(aes(y=..density..), binwidth=0.02, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Basic reproduction number", y = "Density")+
  # geom_vline(xintercept = mean(posterior_samples$R0), col = 2, lwd = 1.2) +
  # geom_vline(xintercept = median(posterior_samples$R0), col = 2, lwd = 1.2, lty = 2) +
  geom_point(aes(x = mean(posterior_samples$R0), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
R0_plot
ggsave("plots/synthetic/uninformative/R0_histogram_uninformative.pdf", plot=R0_plot, device="pdf", 
       width = 6, height = 4)

gamma_plot<-ggplot(posterior_samples, aes(x=gamma)) + 
  geom_histogram(aes(y=..density..), binwidth = 0.0007, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Recovery rate", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$gamma), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
gamma_plot
ggsave("plots/synthetic/uninformative/gamma_histogram_uninformative.pdf", plot=gamma_plot, device="pdf", 
       width = 6, height = 4)


tau_plot<-ggplot(posterior_samples, aes(x=tau)) + 
  geom_histogram(aes(y=..density..),  binwidth = 0.005, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Infection rate", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$tau), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
tau_plot
ggsave("plots/synthetic/uninformative/tau_histogram_uninformative.pdf", plot=tau_plot, 
       device="pdf", width = 6, height = 4)


rho_plot<-ggplot(posterior_samples, aes(x=rho)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Initial proportion of infected", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$rho), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
rho_plot
ggsave("plots/synthetic/uninformative/rho_histogram_uninformative.pdf", plot=rho_plot, 
       device="pdf", width = 6, height = 4)

n_plot<-ggplot(posterior_samples, aes(x=n)) + 
  geom_histogram(aes(y=..density..),  binwidth = 1, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Degree", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$n), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
n_plot
ggsave("plots/synthetic/uninformative/n_histogram_uninformative.pdf", plot=n_plot, 
       device="pdf", width = 6, height = 4)


write.csv(posterior_samples, file="plots/synthetic/uninformative/DSA_posterior_samples_uninformative.csv")



########

params <- c(R0=mean(posterior_samples$R0),
            gamma=mean(posterior_samples$gamma),
            n=mean(posterior_samples$n), 
            rho=mean(posterior_samples$rho))

kT <- length(all_infection_times)
Tmax <- max(all_infection_times) 

effective_population_size <- c()

L <- nrow(posterior_samples)
for (i in 1:L){
  params <- c(R0 = posterior_samples$R0[i],
              gamma = posterior_samples$gamma[i],
              n = posterior_samples$n[i],
              rho = posterior_samples$rho[i])
  xmat=pairwise.ODE(Tmax, params, c(1.0) )
  effective_population_size <- c(effective_population_size, kT/(1- xmat[nrow(xmat),2]))
}

posterior_samples$effective_population_size <- effective_population_size
posterior_samples$I0 <- posterior_samples$effective_population_size * posterior_samples$rho
posterior_samples$total_population_size = posterior_samples$effective_population_size * (1+posterior_samples$rho) 

# hist(effective_population_size)

print(mean(effective_population_size))
print(median(effective_population_size))

I0_plot<-ggplot(posterior_samples, aes(x=I0)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Initial number of infected", y = "Count")+
  geom_point(aes(x = mean(posterior_samples$I0), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
I0_plot
ggsave("plots/synthetic/uninformative/I0_histogram_uninformative.pdf", plot=I0_plot, 
       device="pdf", width = 6, height = 4)


Eff_popsize_plot<-ggplot(posterior_samples, aes(x=effective_population_size)) + 
  geom_histogram(aes(y=..density..), binwidth=10, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Effective population size", y = "Density")+
  # geom_vline(xintercept = mean(posterior_samples$R0), col = 2, lwd = 1.2) +
  # geom_vline(xintercept = median(posterior_samples$R0), col = 2, lwd = 1.2, lty = 2) +
  geom_point(aes(x = mean(posterior_samples$effective_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
Eff_popsize_plot
ggsave("plots/synthetic/uninformative/Effective_popsize_histogram_uninformative.pdf", plot=Eff_popsize_plot, device="pdf", width = 6, height = 4)


Total_popsize_plot<-ggplot(posterior_samples, aes(x=total_population_size)) + 
  geom_histogram(aes(y=..density..), binwidth=10, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Total population size", y = "Density")+
  # geom_vline(xintercept = mean(posterior_samples$R0), col = 2, lwd = 1.2) +
  # geom_vline(xintercept = median(posterior_samples$R0), col = 2, lwd = 1.2, lty = 2) +
  geom_point(aes(x = mean(posterior_samples$total_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
Total_popsize_plot
ggsave("plots/synthetic/uninformative/Total_popsize_histogram_uninformative.pdf", plot=Total_popsize_plot, device="pdf", width = 6, height = 4)


write.csv(posterior_samples, file="plots/synthetic/uninformative/DSA_posterior_samples_uninformative.csv")

