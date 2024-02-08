function [corr_DR_speed_N_d, corr_DR_speed_E_d, corr_DR_lat, corr_DR_long,...
    est_state_new, est_state_err_cov_mat_new]=...
    Integrate_GNSS_DR_soln(est_state_last, est_state_err_cov_mat_last,...
    GNSS_lat_rad, GNSS_long_rad, GNSS_height, GNSS_vel_N, GNSS_vel_E,...
    DR_lat_rad, DR_long_rad, DR_vel_N_d, DR_vel_E_d,...
    prop_time, S_DR, GNSS_pos_sd, GNSS_vel_sd)

% Function to integrate GNSS and DR solutions derived independently
% 
% Inputs : 
% est_state_new : estimated state vector
% est_state_err_cov_mat_last : estimated state error covariance matrix
% 
% GNSS_lat_rad : GNSS derived latitude in radians
% GNSS_long_rad : GNSS derived longitude in radians
% GNSS_heights : GNSS derived height in m
% GNSS_vel_N : GNSS derived velocity in N direction
% GNSS_vel_E : GNSS derived velocity in E direction
% 
% DR_lat_rad : DR latitude solution in radians
% DR_long_rad : DR longitude solution in radians
% DR_vel_N_d : DR velocity (damped) in N direction
% DR_vel_E_d : DR velocity (damped) in E direction
% 
%
% Optional inputs
% 
% prop_time : propogation time (default value in 0.5 s)
% S_DR : Power spectral density (default value in 0.2)
% GNSS_pos_sd : sd of GNSS position solution
% GNSS_vel_sd : sd of GNSS velocity solution
% 
%
% Outputs :
% corr_DR_speed_N_d : Corrected DR speed in N direction
% corr_DR_speed_E_d : Corrected DR speed in E direction
% corr_DR_lat : corrected DR latitude in radians !!!
% corr_DR_long : corrected DR longitude in radians !!!
% est_state_new : estimated state vector for current epoch
% est_state_err_cov_mat_last : estimated state error covariance matrix for currepnt epoch
% 
% 
if ~exist('prop_time','var')
      prop_time = 0.5;
end
if ~exist('S_DR','var')
      S_DR = 0.2;
end
if ~exist('GNSS_pos_sd','var')
      GNSS_pos_sd = 10;
end
if ~exist('GNSS_vel_sd','var')
      GNSS_vel_sd = 0.05;
end
% Kalman filter
% Step 1) 
% define the transition matrix
[R_N, R_E]= Radii_of_curvature(GNSS_lat_rad);
h= GNSS_height;
trans_mat_last= eye(4);
trans_mat_last(3, 1)= prop_time/(R_N+ h );
trans_mat_last(4, 2)= prop_time/((R_E+ h )*cos(GNSS_lat_rad));
% Step 2) 
% System error covariance matrix
sys_noise_cov_mat_last= zeros(4, 4);
sys_noise_cov_mat_last(1,1)= S_DR* prop_time;
sys_noise_cov_mat_last(2,2)= S_DR* prop_time;
sys_noise_cov_mat_last(3,3)= (1/3* S_DR* prop_time^3)/(R_N+ h)^2;
sys_noise_cov_mat_last(4,4)= (1/3* S_DR* prop_time^3) / ((R_E+ h) * cos(GNSS_lat_rad))^2;
sys_noise_cov_mat_last(3, 1)= (1/2 * S_DR* prop_time^2)/(R_N+ h);
sys_noise_cov_mat_last(1, 3)= (1/2 * S_DR* prop_time^2)/(R_N+ h);
sys_noise_cov_mat_last(4, 2)= (1/2 * S_DR* prop_time^2)/((R_E+ h) * cos(GNSS_lat_rad));
sys_noise_cov_mat_last(2, 4)= (1/2 * S_DR* prop_time^2)/((R_E+ h) * cos(GNSS_lat_rad));
% Step 3)
% Propogate state estimates
pred_state_new= trans_mat_last* est_state_last;
% Step 4)
% Propogate error covariance matrix
pred_state_err_cov_mat_new= trans_mat_last* ...
    est_state_err_cov_mat_last* trans_mat_last' + sys_noise_cov_mat_last;
% Step 5)
% Measuerment matrix
meas_mat_new= zeros(4, 4);
meas_mat_new(3:4, 1:2)= -eye(2, 2);
meas_mat_new(1:2, 3:4)= -eye(2, 2);
% Step 6)
% Measurement noise covariance matrix
meas_noise_cov_mat_new= zeros(4, 4);
meas_noise_cov_mat_new(1, 1)= GNSS_pos_sd^2 / (R_N+ h)^2;
meas_noise_cov_mat_new(2, 2)= GNSS_pos_sd^2 / ((R_E+ h)* cos(GNSS_lat_rad))^2;
meas_noise_cov_mat_new(3, 3)= GNSS_vel_sd^2;
meas_noise_cov_mat_new(4, 4)= GNSS_vel_sd^2;
% Step 7)
% Calculate Kalman gain matrix
A= meas_mat_new* pred_state_err_cov_mat_new * meas_mat_new';
B= A+ meas_noise_cov_mat_new;
C= inv(B);
kal_gain_mat_new = pred_state_err_cov_mat_new* meas_mat_new' * C;
% Step 8)
% Forumlate the measurement innovation 
% a) Initialising with GNSS -DR values
pred_meas_inov_new= zeros(4, 1);
pred_meas_inov_new(1) = GNSS_lat_rad- DR_lat_rad;
pred_meas_inov_new(2) = GNSS_long_rad- DR_long_rad;
pred_meas_inov_new(3) = GNSS_vel_N- DR_vel_N_d;
pred_meas_inov_new(4) = GNSS_vel_E- DR_vel_E_d;
% c) subtracting the  measurement innov* state vector
pred_meas_inov_new= pred_meas_inov_new- (meas_mat_new* pred_state_new);
% Step 9) Update state estimates
est_state_new= pred_state_new+ kal_gain_mat_new* pred_meas_inov_new;
% Step 10) Update state error covariance matrix
est_state_err_cov_mat_new= (eye(4)- kal_gain_mat_new* meas_mat_new)* pred_state_err_cov_mat_new;
% Use Kalman filter estimates to correct the DR solution at each epoch
corr_DR_speed_N_d= DR_vel_N_d- est_state_new(1);
corr_DR_speed_E_d= DR_vel_E_d- est_state_new(2);
corr_DR_lat= DR_lat_rad- est_state_new(3);
corr_DR_long= DR_long_rad- est_state_new(4);
end

