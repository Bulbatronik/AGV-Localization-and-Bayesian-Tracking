<h1 align="center">Localization and Bayesian Tracking of an Automatic Guided Vehicle</h1>
*Developed with Esat Ince*.

## Overview

The project required to perform the localization of an automated guided vehicle (AUG) using ultra-wide band signals (UWB). In particular, two different methods of localization should be compared, namely snapshot/static localization techniques and Bayesian Tracking.

## Our approach
### 1. Dataset 
The dataset consisted fo time difference of arrival (TDOA) measurements using 4 UWB tags attached to the vehicle and 6 access points (AP) at known positions. The TDOA measurements were computed wrt. the master AP with id 2. In total, the dataset includes 5 TDOA measurements at sampling rate of 10 Hz. Each TDOA consists of 2000 samples, which corresponds to a 200 seconds of vehicle tracking.              

### 2. Data preprocessing
Given that the measurements come from the real sensors, various outliers and missing values were present. Most importantly, the dataset had no ground truth which could be used to compare the accuracy aof the localization algorithms.

Our first step was to impaint NaN values using the least squares (LS) approach. After removing initial missing values, we applied median filter to identify the identify the outliers, which were replaced using the LS approach. This procedure was repeated three using a large window size of 50, 10 and 3 samples, respectively.

To generate the ground truth, we applied [`Savitzky-Golay` filter](https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter) of order 3 and length 5 to smoothen the data and remove the noise. 

### 3. Methodology
For the localization, we decided to use weighted non-linear least squares (WNLS) for static localization and extended kalman filter (EKF) for Bayesiang tracking. In order to identify the covariance matrix of the measurements, we used the part of the stationary part of the signal, which corresponded to the non-moving vehicle. From this part we extracted the variances of the TDOAs measurements which we used to improve the performance wrt. classic NLS.

EKF is used for Gaussian pdf’s and nonlinear models. We linearize the model around the current location fix and approximate the pdf’s as Gaussians. The most crucial choice is to select the motion model to be used by the EKF. The most common ones are:
- Random walk model
- Nearly constant velocity model
- Nearly constant acceleration model (random jerk)

Random walk model uses a driving process of zero-mean random velocity, which is a suitable model for a pedestrian but not for a car. Nearly Constant Velocity Model and Nearly Constant Acceleration Model are suitable for a car. However, Nearly Constant Velocity model cannot perform well during sharp turns, quick stops or accelerations. We implemented both of the models and showed why Nearly Constant Velocity is performing worse. We also used ellipses to quantify the uncertainty of a-priori and a-posteriori prediction on the AUG's position. 
The final results of the simulation are available in the [Report.pdf](Project/Report.pdf). A quick demo of the tracking is available here:

[Watch the demo video](https://drive.google.com/file/d/1V4K-1m98SW36MdO_wrS5hkyRIjSRFdVu/preview).


## Appendix
The repository also includes a set of lab and homework matirials on usefull algorithms for object localization inside [`Labs`](Labs) folder. The topics include [maximum likelihood (ML) estimation](Labs/1-Static%20Localization/), [lower boiund and iterative methods](Labs/2-Cramer-Rao%20Bound%20and%20iterative%20NLS/), [Kalman filter](Labs/3-Kalman%20filter/) and [Extended Kalman filter](Labs/4-Extended%20Kalman%20Filter/), [Particle filter](Labs/5-Particle%20filter/), as well as an example op [Satellite localization](Labs/6-Satellite%20localization/)  
