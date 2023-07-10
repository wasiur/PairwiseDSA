
setwd("~/OneDrive - The University of Nottingham/Scripts/PairwiseDSA")

#### include packages etc;
options(warn=-1, message =-1)
library(ggplot2); 
library(survival)

library(rstan); 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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
nIter <- 2500 #number of iterations of Hamiltonian Monte Carlo 
nChains <- 8 #number of parallel chains

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


#### Create Stan data

stan_data <- list(N = N, 
                  infection_times = infection_times_sample,
                  M = M,
                  recovery_times = recovery_times_sample)


fit <- stan(
  file = "pairwise_DSA_continuous.stan",
  data = stan_data,
  chains = nChains,
  iter = nIter,
  refresh = 1000
)


stan_summary_fit<-summary(fit)
print(stan_summary_fit)


posterior_samples <- data.frame(a = extract(fit)$a, 
                                b = extract(fit)$b,
                                xi = extract(fit)$xi,
                                rho = extract(fit)$rho,
                                tau = extract(fit)$tau,
                                gamma = extract(fit)$gamma,
                                R0 = extract(fit)$R0,
                                deg = extract(fit)$deg,
                                n=extract(fit)$n)


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
ggsave("plots/synthetic/R0_histogram.pdf", plot=R0_plot, device="pdf", 
       width = 3, height = 6)

gamma_plot<-ggplot(posterior_samples, aes(x=gamma)) + 
  geom_histogram(aes(y=..density..), binwidth = 0.0007, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Recovery rate", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$gamma), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
gamma_plot
ggsave("plots/synthetic/gamma_histogram.pdf", plot=gamma_plot, device="pdf", 
       width = 3, height = 6)


tau_plot<-ggplot(posterior_samples, aes(x=tau)) + 
  geom_histogram(aes(y=..density..),  binwidth = 0.005, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Infection rate", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$tau), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
tau_plot
ggsave("plots/synthetic/tau_histogram.pdf", plot=tau_plot, 
       device="pdf", width = 3, height = 6)


rho_plot<-ggplot(posterior_samples, aes(x=rho)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Initial proportion of infected", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$rho), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
rho_plot
ggsave("plots/synthetic/rho_histogram.pdf", plot=rho_plot, 
       device="pdf", width = 3, height = 6)

n_plot<-ggplot(posterior_samples, aes(x=n)) + 
  geom_histogram(aes(y=..density..),  binwidth = 1, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Degree", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$n), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
n_plot
ggsave("plots/synthetic/n_histogram.pdf", plot=n_plot, 
       device="pdf", width = 3, height = 6)



deg_plot<-ggplot(posterior_samples, aes(x=deg)) + 
  geom_bar(alpha = 0.25, col = 1, fill = 1) +
  labs(x="Degree", y = "Frequency")+
  # geom_point(aes(x = mean(posterior_samples$deg), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
deg_plot
ggsave("plots/synthetic/deg_histogram.pdf", plot=deg_plot, 
       device="pdf", width = 3, height = 6)


write.csv(posterior_samples, file="plots/synthetic/DSA_posterior_samples.csv")



a_plot<-ggplot(posterior_samples, aes(x=a)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="a", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$a), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
a_plot
ggsave("plots/synthetic/a_histogram.pdf", plot=a_plot, 
       device="pdf", width = 3, height = 6)

b_plot<-ggplot(posterior_samples, aes(x=b)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="b", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$b), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
b_plot
ggsave("plots/synthetic/b_histogram.pdf", plot=b_plot, 
       device="pdf", width = 3, height = 6)




xi_plot<-ggplot(posterior_samples, aes(x=xi)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Density of edges", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$xi), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
xi_plot
ggsave("plots/synthetic/xi_histogram.pdf", plot=xi_plot, 
       device="pdf", width = 3, height = 6)

tau_deg_scatter <- ggplot(posterior_samples, aes(x=deg, y=tau)) + geom_point(alpha=0.25, col=1, size=5)+
  labs(x="Degree", y = "Infection rate") + theme_classic()
tau_deg_scatter
ggsave("plots/synthetic/tau_deg_scatter.pdf", plot=tau_deg_scatter, 
       device="pdf", width = 3, height = 6)


ab_scatter <- ggplot(posterior_samples, aes(x=a, y=b)) + geom_point(alpha=0.25, col=1, size=5)+
  labs(x="a", y = "b") + theme_classic()
ab_scatter
ggsave("plots/synthetic/ab_scatter.pdf", plot=ab_scatter, 
       device="pdf", width = 3, height = 6)


## posterior_samples$deg_new = with(posterior_samples, 1+ R0*(tau+gamma)/tau)

## Trace plots
tau_traceplot <- traceplot(fit, "tau")
tau_traceplot
ggsave("plots/synthetic/tau_traceplot.pdf", plot=tau_traceplot, 
       device="pdf", width = 3, height = 6)

gamma_traceplot <- traceplot(fit, "gamma")
gamma_traceplot
ggsave("plots/synthetic/gamma_traceplot.pdf", plot=gamma_traceplot, 
       device="pdf", width = 3, height = 6)

R0_traceplot <- traceplot(fit, "R0")
R0_traceplot
ggsave("plots/synthetic/R0_traceplot.pdf", plot=R0_traceplot, 
       device="pdf", width = 3, height = 6)

xi_traceplot <- traceplot(fit, "xi")
xi_traceplot
ggsave("plots/synthetic/xi_traceplot.pdf", plot=xi_traceplot, 
       device="pdf", width = 3, height = 6)


deg_traceplot <- traceplot(fit, "deg")
deg_traceplot
ggsave("plots/synthetic/deg_traceplot.pdf", plot=deg_traceplot, 
       device="pdf", width = 3, height = 6)

traceplot(fit, "gamma")
traceplot(fit, "rho")
traceplot(fit, "xi")
traceplot(fit, "R0")







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


I0_plot<-ggplot(posterior_samples, aes(x=I0)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Initial number of infected", y = "Count")+
  geom_point(aes(x = mean(posterior_samples$I0), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
I0_plot
ggsave("plots/synthetic/I0_histogram.pdf", plot=I0_plot, 
       device="pdf", width = 3, height = 6)


Eff_popsize_plot<-ggplot(posterior_samples, aes(x=effective_population_size)) + 
  geom_histogram(aes(y=..density..), binwidth=10, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Effective population size", y = "Density")+
  # geom_vline(xintercept = mean(posterior_samples$R0), col = 2, lwd = 1.2) +
  # geom_vline(xintercept = median(posterior_samples$R0), col = 2, lwd = 1.2, lty = 2) +
  geom_point(aes(x = mean(posterior_samples$effective_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
Eff_popsize_plot
ggsave("plots/synthetic/Effective_popsize_histogram.pdf", plot=Eff_popsize_plot, device="pdf", width = 6, height = 4)


Total_popsize_plot<-ggplot(posterior_samples, aes(x=total_population_size)) + 
  geom_histogram(aes(y=..density..), binwidth=10, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Total population size", y = "Density")+
  # geom_vline(xintercept = mean(posterior_samples$R0), col = 2, lwd = 1.2) +
  # geom_vline(xintercept = median(posterior_samples$R0), col = 2, lwd = 1.2, lty = 2) +
  geom_point(aes(x = mean(posterior_samples$total_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
Total_popsize_plot
ggsave("plots/synthetic/Total_popsize_histogram.pdf", plot=Total_popsize_plot, device="pdf", width = 6, height = 4)


write.csv(posterior_samples, file="plots/synthetic/DSA_posterior_samples.csv")


# xmat=pairwise.ODE(Tmax, params, c(1.0) )
# plot(xmat[,1], xmat[,2])

# effective_population_size <- kT/(1- xmat[nrow(xmat),2])
