setwd("~/OneDrive - The University of Nottingham/Scripts/PairwiseDSA")

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



#### Stan parameters 
N <- 5000
nChains <- 4
nIter <- 5000


##### Data 
base_data <- read.csv("data/WSU/Data1.csv")

L <- nrow(base_data)
plot(base_data$reported[1:(nrow(base_data)-1)])



###### small data


Tmax <- 43

# cond1 <- which(base_data$Time > 0)
# cond2 <- which(base_data$Time < Tmax) 
# cond3 <- which(base_data$Time > 0 & base_data$Time < Tmax)
cond1 <- which(1:nrow(base_data) >= 0)
cond2 <- which(1:nrow(base_data) < Tmax)
cond3 <- which(1:nrow(base_data) >= 0 & 1:nrow(base_data) < Tmax)

subset_base_data <- base_data[cond3, , drop = FALSE]




infection_times <- c()

for (i in 1:(L-1)){
  if (base_data$reported[i] > 0){
    # print(i)
    temp <- runif(n = base_data$reported[i], min = i-0.5, max = i+0.5)
    infection_times <- c(infection_times, temp)
  } 
}

hist(infection_times)

infection_times <- infection_times[infection_times < Tmax]

### Take a random sample of infection and recovery times
if (length(infection_times) < N){
  infection_times_sample <- sort(unique(infection_times), decreasing = FALSE)
  N <- length(infection_times_sample)
} else {
  infection_times_sample <- sort(sample(unique(infection_times), size = N, replace = FALSE), decreasing = FALSE)
  #infection_times_sample <- sort(unique(infection_times_sample), decreasing = FALSE)
  N <- length(infection_times_sample)
}


#### Create Stan data

stan_data <- list(N = N, infection_times = infection_times_sample)

# stan_data <- list(N = N, infection_times = infection_times_sample, gamma_0 = gamma_0)



file <- file.path(getwd(),"AH1N1_uninformative_prior.stan")
mod <- cmdstan_model(file)

mod$print()

fit <- mod$sample(
  data = stan_data, 
  seed = 123, 
  chains = nChains, 
  parallel_chains = nChains,
  iter_sampling = nIter,
  refresh = 500 # print update every 500 iters
)

fit$summary()

print(fit$summary())

posterior_samples <- fit$draws(format = "df")

write.csv(posterior_samples, file="plots/AH1N1/small_uninformative/AH1N1_DSA_posterior_samples_uninformative.csv")


R0_plot<-ggplot(posterior_samples, aes(x=R0)) + 
  geom_histogram(aes(y=..density..), binwidth=0.02, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Basic reproduction number", y = "Density")+
  # geom_vline(xintercept = mean(posterior_samples$R0), col = 2, lwd = 1.2) +
  # geom_vline(xintercept = median(posterior_samples$R0), col = 2, lwd = 1.2, lty = 2) +
  geom_point(aes(x = mean(posterior_samples$R0), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
R0_plot
ggsave("plots/AH1N1/small_uninformative/AH1N1_R0_histogram_uninformative.pdf", plot=R0_plot, device="pdf", width = 6, height = 4)

gamma_plot<-ggplot(posterior_samples, aes(x=gamma)) + 
  geom_histogram(aes(y=..density..), binwidth = 0.01, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Recovery rate", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$gamma), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
gamma_plot
ggsave("plots/AH1N1/small_uninformative/AH1N1_gamma_histogram_uninformative.pdf", plot=gamma_plot, device="pdf", width = 6, height = 4)


tau_plot<-ggplot(posterior_samples, aes(x=tau)) + 
  geom_histogram(aes(y=..density..),  binwidth = 0.005, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Infection rate", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$tau), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
tau_plot
ggsave("plots/AH1N1/small_uninformative/AH1N1_tau_histogram_uninformative.pdf", plot=tau_plot, device="pdf", width = 6, height = 4)


rho_plot<-ggplot(posterior_samples, aes(x=rho)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Initial proportion of infected", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$rho), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
rho_plot
ggsave("plots/AH1N1/small_uninformative/AH1N1_rho_histogram_uninformative.pdf", plot=rho_plot, device="pdf", width = 6, height = 4)

n_plot<-ggplot(posterior_samples, aes(x=n)) + 
  geom_histogram(aes(y=..density..),  binwidth = 1, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Degree", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$n), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
n_plot
ggsave("plots/AH1N1/small_uninformative/AH1N1_n_histogram_uninformative.pdf", plot=n_plot, device="pdf", width = 6, height = 4)


params <- c(R0=mean(posterior_samples$R0),
            gamma=mean(posterior_samples$gamma),
            n=mean(posterior_samples$n), 
            rho=mean(posterior_samples$rho))

kT <- length(infection_times)
Tmax <- max(infection_times) 

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
ggsave("plots/AH1N1/small_uninformative/AH1N1_I0_histogram_uninformative.pdf", plot=I0_plot, device="pdf", width = 6, height = 4)


Eff_popsize_plot<-ggplot(posterior_samples, aes(x=effective_population_size)) + 
  geom_histogram(aes(y=..density..),  alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Effective population size", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$effective_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
Eff_popsize_plot
ggsave("plots/AH1N1/small_uninformative/AH1N1_Effective_popsize_histogram_uninformative.pdf", plot=Eff_popsize_plot, device="pdf", width = 6, height = 4)


Total_popsize_plot<-ggplot(posterior_samples, aes(x=total_population_size)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Total population size", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$total_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
Total_popsize_plot
ggsave("plots/AH1N1/small_uninformative/AH1N1_Total_popsize_histogram_uninformative.pdf", plot=Total_popsize_plot, device="pdf", width = 6, height = 4)


write.csv(posterior_samples, file="plots/AH1N1/small_uninformative/AH1N1_DSA_posterior_samples_uninformative.csv")




###### full data
Tmax <- 81

#### Stan parameters 
N <- 1000
nChains <- 1
nIter <- 5000 


##### Data 
base_data <- read.csv("data/WSU/Data1.csv")

L <- nrow(base_data)
plot(base_data$reported[1:(nrow(base_data)-1)])

# cond1 <- which(base_data$Time > 0)
# cond2 <- which(base_data$Time < Tmax) 
# cond3 <- which(base_data$Time > 0 & base_data$Time < Tmax)
cond1 <- which(1:nrow(base_data) >= 0)
cond2 <- which(1:nrow(base_data) < Tmax)
cond3 <- which(1:nrow(base_data) >= 0 & 1:nrow(base_data) < Tmax)

subset_base_data <- base_data[cond3, , drop = FALSE]




infection_times <- c()

for (i in 1:(L-1)){
  if (base_data$reported[i] > 0){
    # print(i)
    temp <- runif(n = base_data$reported[i], min = i-0.5, max = i+0.5)
    infection_times <- c(infection_times, temp)
  } 
}

hist(infection_times)

infection_times <- infection_times[infection_times < Tmax]

### Take a random sample of infection and recovery times
if (length(infection_times) < N){
  infection_times_sample <- sort(unique(infection_times), decreasing = FALSE)
  N <- length(infection_times_sample)
} else {
  infection_times_sample <- sort(sample(unique(infection_times), size = N, replace = FALSE), decreasing = FALSE)
  #infection_times_sample <- sort(unique(infection_times_sample), decreasing = FALSE)
  N <- length(infection_times_sample)
}

hist(infection_times_sample)
#### Create Stan data

stan_data <- list(N = N, infection_times = infection_times_sample)

# stan_data <- list(N = N, infection_times = infection_times_sample, gamma_0 = gamma_0)




file <- file.path(getwd(),"AH1N1_uninformative_prior.stan")
mod <- cmdstan_model(file)

fit <- mod$sample(
  data = stan_data, 
  seed = 123, 
  chains = nChains, 
  parallel_chains = nChains,
  iter_sampling = nIter,
  refresh = 500 # print update every 500 iters
)

fit$summary()

posterior_samples <- fit$draws(format = "df")



R0_plot<-ggplot(posterior_samples, aes(x=R0)) + 
  geom_histogram(aes(y=..density..), binwidth=0.02, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Basic reproduction number", y = "Density")+
  # geom_vline(xintercept = mean(posterior_samples$R0), col = 2, lwd = 1.2) +
  # geom_vline(xintercept = median(posterior_samples$R0), col = 2, lwd = 1.2, lty = 2) +
  geom_point(aes(x = mean(posterior_samples$R0), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
R0_plot
ggsave("plots/AH1N1/full_uninformative/AH1N1_R0_histogram_uninformative.pdf", plot=R0_plot, device="pdf", width = 6, height = 4)

gamma_plot<-ggplot(posterior_samples, aes(x=gamma)) + 
  geom_histogram(aes(y=..density..), binwidth = 0.0007, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Recovery rate", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$gamma), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
gamma_plot
ggsave("plots/AH1N1/full_uninformative/AH1N1_gamma_histogram_uninformative.pdf", plot=gamma_plot, device="pdf", width = 6, height = 4)


tau_plot<-ggplot(posterior_samples, aes(x=tau)) + 
  geom_histogram(aes(y=..density..),  binwidth = 0.001, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Infection rate", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$tau), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
tau_plot
ggsave("plots/AH1N1/full_uninformative/AH1N1_tau_histogram_uninformative.pdf", plot=tau_plot, device="pdf", width = 6, height = 4)


rho_plot<-ggplot(posterior_samples, aes(x=rho)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Initial proportion of infected", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$rho), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
rho_plot
ggsave("plots/AH1N1/full_uninformative/AH1N1_rho_histogram_uninformative.pdf", plot=rho_plot, device="pdf", width = 6, height = 4)

n_plot<-ggplot(posterior_samples, aes(x=n)) + 
  geom_histogram(aes(y=..density..),  binwidth = 1, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Degree", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$n), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
n_plot
ggsave("plots/AH1N1/full_uninformative/AH1N1_n_histogram_uninformative.pdf", plot=n_plot, device="pdf", width = 6, height = 4)



write.csv(posterior_samples, file="plots/AH1N1/full_uninformative/AH1N1_DSA_posterior_samples_uninformative.csv")

########

params <- c(R0=mean(posterior_samples$R0),
            gamma=mean(posterior_samples$gamma),
            n=mean(posterior_samples$n), 
            rho=mean(posterior_samples$rho))

kT <- length(infection_times)
Tmax <- max(infection_times) 

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

hist(effective_population_size)

mean(effective_population_size)
median(effective_population_size)


I0_plot<-ggplot(posterior_samples, aes(x=I0)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, binwidth = 10, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Initial number of infected", y = "Count")+
  geom_point(aes(x = mean(posterior_samples$I0), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic() + 
  xlim(c(0, 1000))
I0_plot
ggsave("plots/AH1N1/full_uninformative/AH1N1_I0_histogram_uninformative.pdf", plot=I0_plot, device="pdf", width = 6, height = 4)


Eff_popsize_plot<-ggplot(posterior_samples, aes(x=effective_population_size)) + 
  geom_histogram(aes(y=..density..),  alpha = 0.25, binwidth = 100, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Effective population size", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$effective_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic() 
Eff_popsize_plot
ggsave("plots/AH1N1/full_uninformative/AH1N1_Effective_popsize_histogram_uninformative.pdf", plot=Eff_popsize_plot, device="pdf", width = 6, height = 4)


Total_popsize_plot<-ggplot(posterior_samples, aes(x=total_population_size)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Total population size", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$total_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic() + 
  xlim(c(0, 40000))
Total_popsize_plot
ggsave("plots/AH1N1/full_uninformative/AH1N1_Total_popsize_histogram_uninformative.pdf", plot=Total_popsize_plot, device="pdf", width = 6, height = 4)


write.csv(posterior_samples, file="plots/AH1N1/full_uninformative/AH1N1_DSA_posterior_samples_uninformative.csv")

abc <- read.csv(file="plots/AH1N1/full_uninformative/AH1N1_DSA_posterior_samples_uninformative.csv")

Eff_popsize_plot<-ggplot(abc, aes(x=effective_population_size)) + 
  geom_histogram(aes(y=..density..),  alpha = 0.25, binwidth = 100, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Effective population size", y = "Density")+
  geom_point(aes(x = mean(abc$effective_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic() 
Eff_popsize_plot
ggsave("plots/AH1N1/full_uninformative/AH1N1_Effective_popsize_histogram_uninformative.pdf", plot=Eff_popsize_plot, device="pdf", width = 6, height = 4)


