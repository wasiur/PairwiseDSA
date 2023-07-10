setwd("~/OneDrive - The University of Nottingham/Scripts/PairwiseDSA")

#### include packages etc;
options(warn=-1, message =-1)
library(ggplot2); 
library(survival)

library(rstan); 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("tudColours.R")

library(deSolve)


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


## Stan configurations
N <- 5000 #sample size for infection data
M <- 5000 #sample size for recovery data 
nIter <- 5000 #number of iterations of Hamiltonian Monte Carlo 
nChains <- 4 #number of parallel chains



##### Data 
base_data <- read.csv("data/COVID19/case_time_series.csv", 
                      colClasses=c("character", "Date", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric" ))


plot(base_data$Daily.Confirmed[1:(nrow(base_data))])

start_date <- as.Date("2021-02-15")
end_date <- as.Date("2021-06-30")

cond1 <- which(base_data$Date_YMD > start_date)
cond2 <- which(base_data$Date_YMD < end_date) 
cond3 <- which(base_data$Date_YMD > start_date & base_data$Date_YMD < end_date)
# cond1 <- which(1:nrow(base_data) >= 600)
# cond2 <- which(1:nrow(base_data) < 400)
# cond3 <- which(1:nrow(base_data) >= 400 & 1:nrow(base_data) < 600)
df <- data.frame(x = base_data$Date_YMD,
  # x = 1:nrow(base_data),
  y = base_data$Daily.Confirmed[1:(nrow(base_data))])

options(scipen=100)
covid_cases <- ggplot() + 
  theme_classic() +
  xlab("Date")+
  ylab("Daily new cases")+
  geom_point(data = df[cond1,],
            mapping = aes(x = x, y = y), lwd = 1.0) +
  geom_point(data = df[cond2,],
            mapping = aes(x = x, y = y), lwd = 1.0) + 
  geom_point(data = df[cond3,],
            mapping = aes(x = x, y = y), col = 2,  lwd = 1.0) +
  geom_rect(aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
            alpha = 0.25, fill = 2) 
  # ggtitle(label = "A cool title", subtitle = "with subtitle")

covid_cases
ggsave("plots/COVID19/covid_cases.pdf", plot=covid_cases, device="pdf", width = 6.0, height = 3)


subset_base_data <- base_data[which(difftime(base_data$Date_YMD, start_date) >0), , drop = FALSE]
subset_base_data<- subset_base_data[which(difftime(subset_base_data$Date_YMD, end_date) <0), ,drop = FALSE]

subset_base_data

L <- nrow(subset_base_data)
infection_times <- c()

for (i in 1:(L-1)){
  if (subset_base_data$Daily.Confirmed[i] > 0){
    # print(i)
    temp <- runif(n = subset_base_data$Daily.Confirmed[i], min = i-0.5, max = i+0.5)
    infection_times <- c(infection_times, temp)
  } 
}

plot(subset_base_data$Daily.Confirmed[1:(nrow(subset_base_data)-1)])
hist(infection_times)

Tmax <- 120
infection_times <- infection_times[infection_times < Tmax]


recovery_times <- c()
for (i in 1:(L-1)){
  if (subset_base_data$Daily.Recovered[i] > 0){
    # print(i)
    temp <- runif(n = subset_base_data$Daily.Recovered[i], min = i-0.5, max = i+0.5)
    recovery_times <- c(recovery_times, temp)
  } 
}
hist(recovery_times)

### Take a random sample of infection and recovery times
if (length(infection_times) < N){
  infection_times_sample <- sort(unique(infection_times), decreasing = FALSE)
  N <- length(infection_times_sample)
} else {
  infection_times_sample <- sort(sample(unique(infection_times), size = N, replace = FALSE), decreasing = FALSE)
  #infection_times_sample <- sort(unique(infection_times_sample), decreasing = FALSE)
  N <- length(infection_times_sample)
}

if (length(recovery_times) < M){
  recovery_times_sample <- sort(unique(recovery_times), decreasing = FALSE)
  M <- length(recovery_times_sample)
} else {
  recovery_times_sample <- sort(sample(unique(recovery_times), size = M, replace = FALSE), decreasing = FALSE)
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



# posterior_samples <- read.csv("plots/COVID19/COVID19_DSA_posterior_samples.csv")
# summary_table <- summary(posterior_samples)
# write.csv(summary_table, file = "plots/COVID19/COVID19_posterior_summary.csv")



R0_plot<-ggplot(posterior_samples, aes(x=R0)) + 
  geom_histogram(aes(y=..density..), binwidth=0.02, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Basic reproduction number", y = "Density")+
  # geom_vline(xintercept = mean(posterior_samples$R0), col = 2, lwd = 1.2) +
  # geom_vline(xintercept = median(posterior_samples$R0), col = 2, lwd = 1.2, lty = 2) +
  geom_point(aes(x = mean(posterior_samples$R0), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
R0_plot
ggsave("plots/COVID19/COVID19_R0_histogram.pdf", plot=R0_plot, device="pdf", width = 6, height = 4)

gamma_plot<-ggplot(posterior_samples, aes(x=gamma)) + 
  geom_histogram(aes(y=..density..), binwidth = 0.0007, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Recovery rate", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$gamma), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
gamma_plot
ggsave("plots/COVID19/COVID19_gamma_histogram.pdf", plot=gamma_plot, device="pdf", width = 6, height = 4)


tau_plot<-ggplot(posterior_samples, aes(x=tau)) + 
  geom_histogram(aes(y=..density..),  binwidth = 0.005, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Infection rate", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$tau), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
tau_plot
ggsave("plots/COVID19/COVID19_tau_histogram.pdf", plot=tau_plot, device="pdf", width = 6, height = 4)


rho_plot<-ggplot(posterior_samples, aes(x=rho)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Initial proportion of infected", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$rho), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
rho_plot
ggsave("plots/COVID19/COVID19_rho_histogram.pdf", plot=rho_plot, device="pdf", width = 6, height = 4)

n_plot<-ggplot(posterior_samples, aes(x=n)) + 
  geom_histogram(aes(y=..density..),  binwidth = 1, alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Degree", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$n), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
n_plot
ggsave("plots/COVID19/COVID19_n_histogram.pdf", plot=n_plot, device="pdf", width = 6, height = 4)



deg_plot<-ggplot(posterior_samples, aes(x=deg)) + 
  geom_bar(alpha = 0.25, col = 1, fill = 1) +
  labs(x="Degree", y = "Frequency")+
  # geom_point(aes(x = mean(posterior_samples$deg), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
deg_plot
ggsave("plots/COVID19/COVID19_deg_histogram.pdf", plot=deg_plot, device="pdf", width = 6, height = 4)


write.csv(posterior_samples, file="plots/COVID19/COVID19_DSA_posterior_samples.csv")



a_plot<-ggplot(posterior_samples, aes(x=a)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="a", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$a), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
a_plot
ggsave("plots/COVID19/COVID19_a_histogram.pdf", plot=a_plot, device="pdf", width = 6, height = 4)

b_plot<-ggplot(posterior_samples, aes(x=b)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="b", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$b), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
b_plot
ggsave("plots/COVID19/COVID19_b_histogram.pdf", plot=b_plot, device="pdf", width = 6, height = 4)




xi_plot<-ggplot(posterior_samples, aes(x=xi)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Density of edges", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$xi), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
xi_plot
ggsave("plots/COVID19/COVID19_xi_histogram.pdf", plot=xi_plot, device="pdf", width = 6, height = 4)

tau_deg_scatter <- ggplot(posterior_samples, aes(x=deg, y=tau)) + geom_point(alpha=0.25, col=1, size=5)+
  labs(x="Degree", y = "Infection rate") + theme_classic()
tau_deg_scatter
ggsave("plots/COVID19/COVID19_tau_deg_scatter.pdf", plot=tau_deg_scatter, device="pdf", width = 6, height = 4)


ab_scatter <- ggplot(posterior_samples, aes(x=a, y=b)) + geom_point(alpha=0.25, col=1, size=5)+
  labs(x="a", y = "b") + theme_classic()
ab_scatter
ggsave("plots/COVID19/COVID19_ab_scatter.pdf", plot=ab_scatter, device="pdf", width = 6, height = 4)


## posterior_samples$deg_new = with(posterior_samples, 1+ R0*(tau+gamma)/tau)

## Trace plots
tau_traceplot <- traceplot(fit, "tau")
tau_traceplot
ggsave("plots/COVID19/COVID19_tau_traceplot.pdf", plot=tau_traceplot, device="pdf", width = 6, height = 4)

gamma_traceplot <- traceplot(fit, "gamma")
gamma_traceplot
ggsave("plots/COVID19/COVID19_gamma_traceplot.pdf", plot=gamma_traceplot, device="pdf", width = 6, height = 4)

R0_traceplot <- traceplot(fit, "R0")
R0_traceplot
ggsave("plots/COVID19/COVID19_R0_traceplot.pdf", plot=R0_traceplot, device="pdf", width = 6, height = 4)

xi_traceplot <- traceplot(fit, "xi")
xi_traceplot
ggsave("plots/COVID19/COVID19_xi_traceplot.pdf", plot=xi_traceplot, device="pdf", width = 6, height = 4)


deg_traceplot <- traceplot(fit, "deg")
deg_traceplot
ggsave("plots/COVID19/COVID19_deg_traceplot.pdf", plot=deg_traceplot, device="pdf", width = 6, height = 4)

traceplot(fit, "gamma")
traceplot(fit, "rho")
traceplot(fit, "xi")
traceplot(fit, "R0")







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

# hist(effective_population_size)


I0_plot<-ggplot(posterior_samples, aes(x=I0)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.75)+
  labs(x="Initial number of infected", y = "Count")+
  geom_point(aes(x = mean(posterior_samples$I0), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
I0_plot
ggsave("plots/COVID19/COVID19_I0_histogram.pdf", plot=I0_plot, device="pdf", width = 6, height = 4)


Eff_popsize_plot<-ggplot(posterior_samples, aes(x=effective_population_size)) + 
  geom_histogram(aes(y=..density..),  alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Effective population size", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$effective_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
Eff_popsize_plot
ggsave("plots/COVID19/COVID19_Effective_popsize_histogram.pdf", plot=Eff_popsize_plot, device="pdf", width = 6, height = 4)


Total_popsize_plot<-ggplot(posterior_samples, aes(x=total_population_size)) + 
  geom_histogram(aes(y=..density..), alpha = 0.25, col = 1, fill = 1) +
  geom_density(alpha=.8, lwd = 1.2, adjust = 1.5)+
  labs(x="Total population size", y = "Density")+
  geom_point(aes(x = mean(posterior_samples$total_population_size), y = 0), pch = 17, size = 10, col = 2) +
  theme_classic()
Total_popsize_plot
ggsave("plots/COVID19/COVID19_Total_popsize_histogram.pdf", plot=Total_popsize_plot, device="pdf", width = 6, height = 4)


write.csv(posterior_samples, file="plots/COVID19/COVID19_DSA_posterior_samples.csv")


# xmat=pairwise.ODE(Tmax, params, c(1.0) )
# plot(xmat[,1], xmat[,2])

# effective_population_size <- kT/(1- xmat[nrow(xmat),2])

print(stan_summary_fit)


#####
plot(subset_base_data$Date_YMD, subset_base_data$Total.Confirmed)
posterior_samples <- read.csv("plots/COVID19/COVID19_DSA_posterior_samples.csv")

nDays <- nrow(subset_base_data)
nSim <- nrow(posterior_samples)
fitted_cum_curves <- array(data=NA, dim = c(nDays,nSim))
# fitted_curves <- subset_base_data
fitted_daily <- array(data=NA, dim = c(nDays,nSim))


for (i in 1:nSim) {
  params <- c(R0 = posterior_samples$R0[i],
              gamma = posterior_samples$gamma[i],
              n = posterior_samples$n[i],
              rho = posterior_samples$rho[i])
  times <- seq(1, nDays, by = 1.0)
  out <- ode(func=S_dot, y=c(S=1.0), times=times, parms=params) 
  n_eff = posterior_samples$effective_population_size[i]
  column_name <- paste0("simulated_cum_infected_", i)
  fitted_cum_curves[, i] <- n_eff*(1-out[,2]) + (subset_base_data$Total.Confirmed[1]-subset_base_data$Daily.Confirmed[1])
  # fitted_curves[column_name] <- n_eff*(1-out[,2]) + (fitted_curves$Total.Confirmed[1]-fitted_curves$Daily.Confirmed[1])
  # column_name <- paste0("simulated_daily_", i)
  fitted_daily[, i] <- c(subset_base_data$Daily.Confirmed[1] ,diff(fitted_cum_curves[, i]))
}

mean_cum_curve <- rowMeans(fitted_cum_curves)
std_cum_curve <- apply(fitted_cum_curves, 1, sd)
low_cum_curve <- mean_cum_curve - 1.96*std_cum_curve
high_cum_curve <- mean_cum_curve + 1.96*std_cum_curve

mean_daily_curve <- rowMeans(fitted_daily)
std_daily_curve <- apply(fitted_daily, 1, sd)
low_daily_curve <- mean_daily_curve - 1.96*std_daily_curve
high_daily_curve <- mean_daily_curve + 1.96*std_daily_curve


covid_fitted_cumulative <- ggplot(data = subset_base_data, aes(x=subset_base_data$Date_YMD, y=subset_base_data$Total.Confirmed), col = 2) +
  theme_classic()+
  xlab("Date")+
  ylab("Cumulative cases")+
  geom_line(aes(x=subset_base_data$Date_YMD, y=subset_base_data$Total.Confirmed), col = 1) + 
  geom_line(aes(x=subset_base_data$Date_YMD, y=mean_cum_curve), col = 2) + 
  geom_ribbon(aes(ymin=low_cum_curve, ymax=high_cum_curve), col = 4, alpha = 0.25)

covid_fitted_cumulative
ggsave(covid_fitted_cumulative, filename = "plots/COVID19/COVID19_fitted_cumulative_cases.pdf", device="pdf", width = 6, height = 4)

covid_fitted_daily <- ggplot(data = subset_base_data, aes(x=subset_base_data$Date_YMD, y=subset_base_data$Daily.Confirmed), col = 2) +
  theme_classic()+
  xlab("Date")+
  ylab("Daily new cases")+
  geom_line(aes(x=subset_base_data$Date_YMD, y=subset_base_data$Daily.Confirmed), col = 1) + 
  geom_line(aes(x=subset_base_data$Date_YMD, y=mean_daily_curve), col = 2) + 
  geom_ribbon(aes(ymin=low_daily_curve, ymax=high_daily_curve), col = 4, alpha = 0.25)

covid_fitted_daily
ggsave(covid_fitted_daily, filename = "plots/COVID19/COVID19_fitted_daily_cases.pdf", device="pdf", width = 6, height = 4)


write.csv(fitted_cum_curves, "plots/COVID19/COVID19_fitted_cum_cases.csv")
write.csv(fitted_daily, "plots/COVID19/COVID19_fitted_daily_cases.csv")



library(zoo)
ma_length <- 14
ma_cum_curve <- rollmean(subset_base_data$Total.Confirmed, ma_length)
ma_variance_factor <- sd(ma_cum_curve-subset_base_data$Total.Confirmed[ma_length:nrow(subset_base_data)])


mean_cum_curve <- rowMeans(fitted_cum_curves)
std_cum_curve <- apply(fitted_cum_curves, 1, sd)
low_cum_curve <- mean_cum_curve - ma_variance_factor*1.96*std_cum_curve
high_cum_curve <- mean_cum_curve + ma_variance_factor*1.96*std_cum_curve


ma_daily_curve <- rollmean(subset_base_data$Daily.Confirmed, ma_length)
ma_variance_factor <- sd(ma_daily_curve-subset_base_data$Daily.Confirmed[ma_length:nrow(subset_base_data)])

mean_daily_curve <- rowMeans(fitted_daily)
std_daily_curve <- apply(fitted_daily, 1, sd)
low_daily_curve <- mean_daily_curve - ma_variance_factor*1.96*std_daily_curve
high_daily_curve <- mean_daily_curve + ma_variance_factor*1.96*std_daily_curve


COVID19_fitted_cumulative <- ggplot(data = subset_base_data, aes(x=subset_base_data$Date_YMD, y=subset_base_data$Total.Confirmed), col = 2) +
  theme_classic()+
  xlab("Time")+
  ylab("Cumulative cases")+
  geom_line(aes(x=subset_base_data$Date_YMD, y=subset_base_data$Total.Confirmed), col = 1) + 
  geom_line(aes(x=subset_base_data$Date_YMD, y=mean_cum_curve), col = 2) + 
  geom_ribbon(aes(ymin=low_cum_curve, ymax=high_cum_curve), col = 4, alpha = 0.25)

COVID19_fitted_cumulative
ggsave(fmd_fitted_cumulative, filename = "plots/COVID19/COVID19_fitted_cumulative_cases_variance_adjusted.pdf", device="pdf", width = 6, height = 4)


COVID19_fitted_daily <- ggplot(data = subset_base_data, aes(x=subset_base_data$Date_YMD, y=subset_base_data$Daily.Confirmed), col = 2) +
  theme_classic()+
  xlab("Date")+
  ylab("Daily new cases")+
  geom_line(aes(x=subset_base_data$Date_YMD, y=subset_base_data$Daily.Confirmed), col = 1) + 
  geom_line(aes(x=subset_base_data$Date_YMD, y=mean_daily_curve), col = 2) + 
  geom_ribbon(aes(ymin=low_daily_curve, ymax=high_daily_curve), col = 4, alpha = 0.25)

COVID19_fitted_daily
ggsave(fmd_fitted_daily, filename = "plots/COVID19/COVID19_fitted_daily_cases_variance_adjusted.pdf", device="pdf", width = 6, height = 4)



