function [est_state_new, est_error_cov_mat_new]= Kalman_filter_basic(est_state_last,...
    est_error_cov_mat_last, GNSS_pos, tau_s)

% function to implement Basic Kalman filter
% 
% Inputs :
% est_state_last : estimated state vector from previous epoch (6, 1)
% est_error_cov_mat_last : estimated error covariance matrix from previous epoch (6, 6)
% tau_s : delta t between previous and current epoch (1, 1)
% 
% Returns:
% est_state_new : estimated state vector for current epoch (6, 1)
% est_error_cov_mat_new : estimated error covariance matrix for current epoch (6, 6)


% transition matrix as a function of delta_t 
trans_mat_last = zeros(6, 6);
trans_mat_last(1:3, 1:3) = eye(3);
trans_mat_last(4:6, 1:3) = zeros(3, 3);
trans_mat_last(1:3, 4:6) = eye(3)* tau_s;
trans_mat_last(4:6, 4:6) = eye(3);


% system noise covariave matrix
noise_cov_mat_last= zeros(6 ,6);
noise_cov_mat_last(1:3, 1:3)= (1/3)*5*(tau_s^3)*eye(3);
noise_cov_mat_last(1:3, 4:6)= (1/2)*5*(tau_s^2)*eye(3);
noise_cov_mat_last(4:6, 1:3)= (1/2)*5*(tau_s^2)*eye(3);
noise_cov_mat_last(4:6, 4:6)= 5*tau_s*eye(3);


% propogate the state vector
pred_state_new= trans_mat_last* est_state_last;





% propogate the error covariance matrix
pred_error_cov_mat_new= trans_mat_last* est_error_cov_mat_last* trans_mat_last' + noise_cov_mat_last;






%  Formulate the measurement matrix
meas_mat_new= zeros(3, 6);
meas_mat_new(1:3, 1:3) = eye(3);



% Compute measurement noise covariance matrix
meas_noise_cov_mat_new= eye(3)* 2.5^2;



% Compute the Kalman gain matrix
A= meas_mat_new* pred_error_cov_mat_new * meas_mat_new';
B= A+ meas_noise_cov_mat_new;
C= inv(B);
kal_gain_mat_new = pred_error_cov_mat_new* meas_mat_new' * C;





% Formulate measurment innovation vector
pred_meas_inov_new= GNSS_pos' - pred_state_new(1:3, 1);



% Update state vector
est_state_new= pred_state_new+ kal_gain_mat_new* pred_meas_inov_new;


% Update error covariance matrix
A= kal_gain_mat_new* meas_mat_new ;
B= eye(size(A, 1));

est_error_cov_mat_new = (B- A)* pred_error_cov_mat_new ;





