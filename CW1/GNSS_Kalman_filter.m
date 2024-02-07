function [est_state_new, est_error_cov_mat_new]= ...
    GNSS_Kalman_filter(est_state_last, est_error_cov_mat_last,...
    sat_index, sat_count, simulation_time, tau_s,...
    meas_range_asati, meas_range_rate_asati)

% Function to implement GNSS based Kalman filter
% Incorporates the pseudo range and pseudo range rate measurements from all satellites
% 
% Inputs :
% est_state_last : estimated state vector of k-1 epoch
% est_error_cov_mat_last : estimated error covariance matrix of k-1 epoch
% sat_index : list of index of all satellites from which measurements are available
% sat_count : lenght of sat_index
% simulation_time : time in seconds (required to get correct position and velocity of all satellites)
% tau_s : propogation interval in seconds (required to compute the transition and system noise covariance matrix
% meas_range_asati : pseudo range measurements from all satellites
% meas_range_rate_asati : pseudo range rate measurements from all satellites
% 
% 
% Returns:
% est_state_new : estimated state vecor for epoch k
% est_error_cov_mat_new : estimated error covariance matrix at epoch k


if size(meas_range_asati, 1) > size(meas_range_asati, 2)
    % no. of rows > no. of col
    % col vector
    % NOT compatible with the rest of the code
    meas_range_asati = meas_range_asati';
end

if size(meas_range_rate_asati, 1) > size(meas_range_rate_asati, 2)
    % no. of rows > no. of col
    % col vector
    % NOT compatible with the rest of the code
    meas_range_rate_asati = meas_range_rate_asati';
end


% if ~exist('prop_time','var')
%       prop_time = 0.5;
% end





% Only some constatns from Define_Constants.m are needed
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
Omega_ie = Skew_symmetric([0,0,omega_ie]);
c = 299792458; % Speed of light in m/s





% Kalman filter Step 1)
% Transisition matrix last

trans_mat_last = eye(8, 8);


trans_mat_last(1:3, 4:6)= eye(3)* tau_s;
trans_mat_last(7, 8)= tau_s;




% Kalman filter Step 2) 
% Compute system noise covariance matrix


acc_pow_spec_dens= 5; % can be added as function argument
clock_phase_dens= 0.01; % can be added as function argument
clock_freq_dens= 0.04; % can be added as function argument


noise_cov_mat_last= zeros(8 , 8);


noise_cov_mat_last(1:3, 1:3)= (1/3)*acc_pow_spec_dens*(tau_s^3)*eye(3);

noise_cov_mat_last(1:3, 4:6)= (1/2)*acc_pow_spec_dens*(tau_s^2)*eye(3);

noise_cov_mat_last(4:6, 1:3)= (1/2)*acc_pow_spec_dens*(tau_s^2)*eye(3);

noise_cov_mat_last(4:6, 4:6)= acc_pow_spec_dens*tau_s*eye(3);

noise_cov_mat_last(7, 7) = clock_phase_dens*tau_s + (1/3)*clock_freq_dens*(tau_s^3) ;
noise_cov_mat_last(7, 8) = (1/2)*clock_freq_dens*(tau_s^2);
noise_cov_mat_last(8, 7) = (1/2)*clock_freq_dens*(tau_s^2);
noise_cov_mat_last(8, 8) = clock_freq_dens*tau_s;




% Step 3)
% propogate state estimates pred_state_new using transisiton matrix

pred_state_new= trans_mat_last* est_state_last;



% Step 4)
% pred_error_cov_mat_new
% propogate the error covariance matrix

pred_error_cov_mat_new= trans_mat_last* est_error_cov_mat_last* trans_mat_last' + noise_cov_mat_last;



% Getting position and velocity of all satellites

pred_r_esati= zeros(sat_count, 3);
pred_v_esati= zeros(sat_count, 3);


for index = 1: numel(sat_index)
    [pred_r_esati(index, :) , pred_v_esati(index, :)]= ...
    Satellite_position_and_velocity(simulation_time, sat_index(index));
end


% Compute psudo ranges and line-of-sight unit vector from antenna to all satellites

pred_range_asati= zeros(sat_count, 1);
pred_u_asati= zeros(sat_count, 3);

for index=1:sat_count
   [pred_range_asati(index, 1), pred_u_asati(index, :)]= ...
       line_of_sight_vector(pred_r_esati(index, :)', pred_state_new(1:3, 1));

end



% Compute range rates from approx antenna position to satellite

pred_range_rate_asati= zeros(sat_count, 1);

for index=1:sat_count

    C_I_e= eye(3); % (3, 3)
    C_I_e(1, 2)= omega_ie*pred_range_asati(index, 1)/c;
    C_I_e(2, 1)= -omega_ie*pred_range_asati(index, 1)/c;
    
    A= pred_v_esati(index, :)'+ Omega_ie*pred_r_esati(index, :)'; % (3, 1)

    B= pred_state_new(4:6, 1)+ Omega_ie* pred_state_new(1:3, 1); % (3, 1)

    pred_range_rate_asati(index, 1)= pred_u_asati(index, :)*((C_I_e* A)- B);
end



% Step 5) Compute measurement matrix

meas_mat_new = zeros(sat_count*2, 8);

meas_mat_new(1:sat_count, 1:3)= -pred_u_asati;
meas_mat_new(1:sat_count, 7)= 1;

meas_mat_new(sat_count+1: sat_count*2, 4:6)= -pred_u_asati;
meas_mat_new(sat_count+1: sat_count*2, 8)= 1;



% Step 6)
%  Compute measurement noise covariance matrix

pseudo_range_sd= 10;
pseudo_range_rate_sd= 0.05;


meas_noise_cov_mat_new= eye(sat_count*2, sat_count*2);

meas_noise_cov_mat_new(1: sat_count, 1: sat_count)= eye(sat_count, sat_count)* pseudo_range_sd^2;

meas_noise_cov_mat_new(sat_count+1: 2*sat_count, sat_count+1: 2*sat_count)= ...
    eye(sat_count, sat_count)* pseudo_range_rate_sd^2;


% Step 7)
% Compute kalman gain matrix
%  Kalman_gain_mat_new

A= meas_mat_new* pred_error_cov_mat_new * meas_mat_new';

B= A+ meas_noise_cov_mat_new;

C= pinv(B);



kal_gain_mat_new = pred_error_cov_mat_new* meas_mat_new' * C;



% Step 8)
% Predicted measurement innovation vector pred_meas_inov_new
pred_meas_inov_new= zeros(sat_count*2, 1);

pred_meas_inov_new(1: sat_count, 1)= meas_range_asati'- pred_range_asati- pred_state_new(7, 1);
pred_meas_inov_new(sat_count+1: 2*sat_count, 1)= ...
    meas_range_rate_asati'- pred_range_rate_asati- pred_state_new(8, 1);


% Step 9)
%  Update state vector

est_state_new= pred_state_new+ kal_gain_mat_new* pred_meas_inov_new;


% Step 10)
%  Update error covariance matrix

A= kal_gain_mat_new* meas_mat_new ;
B= eye(size(A, 1));

est_error_cov_mat_new = (B- A)* pred_error_cov_mat_new;


end
