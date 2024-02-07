function [x_est,P_matrix] = Initialise_GNSS_KF_CW
%Initialise_GNSS_KF - Initializes the GNSS EKF state estimates and error
%covariance matrix for Workshop 2
%
% This function created 30/11/2016 by Paul Groves
%
% Outputs:
%   x_est                 Kalman filter estimates:
%     Rows 1-3            estimated ECEF user position (m)
%     Rows 4-6            estimated ECEF user velocity (m/s)
%     Row 7               estimated receiver clock offset (m) 
%     Row 8               estimated receiver clock drift (m/s)
%   P_matrix              state estimation error covariance matrix

% Copyright 2016, Paul Groves
% License: BSD; see license.txt for details

% Begins

% Initialise state estimates 
x_est = [  3.977851131656157e+06; -1.118086981210384e+04; 4.969033742600406e+06 ;...% positon
    0.006799278282820; 0.047183610241202;  -0.017330817394829;...% velocity
    1.000880095199524e+04; 1.000123511027067e+02]; % clock offset , clock drift

% Initialise error covariance matrix
P_matrix =  zeros(8);
P_matrix(1,1) = 10^2;
P_matrix(2,2) = 10^2;
P_matrix(3,3) = 10^2;
P_matrix(4,4) = 0.05^2;
P_matrix(5,5) = 0.05^2;
P_matrix(6,6) = 0.05^2;
P_matrix(7,7) = 100000^2;
P_matrix(8,8) = 200^2;

% Ends