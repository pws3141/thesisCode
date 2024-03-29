
# Thesis Code

In this git repository, we have included some R code used thesis "Forecasting
Complex System Using Stochastic Models for Low Dimensional Approximations".
This code is split into folders mirroring the layout of the thesis.

# Part I: Theory

This section of the thesis introduces the three main themes: dimenion reduction,
stochastic modelling, and forecasting.

## Ch1: Dimension Reduction

In this section we discussed PCA and ICA. Within ICA, the fastICA method was
critiqued, and the (novel) clusterICA method was introduced.

### fastICA

This includes the code to reproduce the figures in the fastICA section (1.4.2). This
section also forms a paper "On the estimation of entropy in the fastICA method",
which has been submitted and is currently in the review stage.

### clusterICA

This section includes the clusterICA code.

## Ch2: Modelling using Stochastic Processes

This code avaliable here intoduced in this section is taken from the parameter
fitting procedures (Section 2.3.2), and estimating the seasonal effects (Section
2.3.3).

### Ornstein-Uhlenbeck parameter estimation

Here we include the maximum likelihood estimation for sparse observations, alongside
the method of moments.

### Block-average Ornstein-Uhlenbeck parameter estimation

Here we include the maximum likelihood estimation method used to obtain drift
and diffusion parameter estimations from a block-average realisation.

### Seasonal effects

In this section we include the code to obtain a quadratic spline given some seasonal
means.

## Ch3: Forecasting

Here we include the code that produces an example of point forecast and prediction intervals for
a realisation from Ornstein-Uhlenbeck process. We also include the code used to calculate
one-step-ahead point forecasts and prediction intervals for block-average
Ornstein-Uhlenbeck processes.

# Part II: Application

## Ch4: HadCM3: mean sea-level air pressure

## Ch5: HadCM3: mean sea-level air pressure with wind velocity

Here we include code for reconstructing the climate simulator output.

### Plotting techniques

This code creates the plots show in Section 5.5, including the arrows representing
the wind velocity.
