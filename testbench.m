clc; clear all; close all; format compact;
%
%
% Matlab version of the algorithm by Wolf et al. for estimating the dominant Lyapunov exponent from a 1-D time series.
%
% Physica 16D (1985) 285-317 "Determining Lyapunov Exponents from a Time Series"
% Alan Wolf, Jack B. Swift, Harry L. Swinney, and John A. Vastano
% 
% Appendix B of the Physica D article contains Fortran code for a concise, but highly inefficient version of the algorithm.
% I have been distributing a Fortran and C version of the efficient version of the algorithm since the 1980's.
% The efficient version of the code was converted to Matlab by Taehyeun Park, The Cooper Union, EE'15 in September, 2014.
%
% Detailed instructions for the use of this code will be posted at Matlab Central's File Exchange.
% 
% This file, testbench.m, takes 16,384 points of Lorenz attractor data (dt = 0.01 seconds) and estimates the dominant 
% Lyapunov exponent by first calling "basgen" (preparing a database for quickly finding nearest neighbors in reconstructed
% phase space) and then calling "fet" to estimate the Lyapunov exponent.
% The program shows both graphical output as orbital divergence is monitored and a text file showing the running exponent estimate.
%
%

fname = 'Data.lor';

datcnt = 16384;
tau = 10;
ndim = 3;
ires = 10;
maxbox = 6000;

db = basgen(fname, tau, ndim, ires, datcnt, maxbox);

dt = .01;
evolve = 20;
dismin = 0.001;
dismax = 0.3;
thmax = 30;

[out, SUM] = fet(db, dt, evolve, dismin, dismax, thmax);

makeplot(db, out, evolve, 'NorthWest')