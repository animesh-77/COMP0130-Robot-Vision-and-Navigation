function [est_state] = Single_Epoch_position(time, ...
    pred_pos_e_ea, ~, sat_index, obs_prange_a_sati)
% 
% 
% Function to get single epoch position solution with OR without any priors
% Iterates until position solution does not improve by more than 10cm in x y z axis in ECEF coordinate system
% 
% 
% 
% Inputs :
% time : time in seconds (1, 1)
% pred_pos_e_ea : prior position of the antenna. If not available then coordinates for centre of the Earth is passed (3, 1)  
% ~ : velocity vector not required (3, 1)
% sat_index : index number of satellites (n,) where n is the number of satellites
% obs_prange_a_sati : observed pseudo range measurements from all satellites (n, 1)
% 
% 
% 
% Returns
% est_state : position of the antenna in ECEF coordinate system and clock
% offset
Define_Constants
pred_clock_offset= 1;
sat_count= numel(sat_index);
% total number of satellites
est_pos_e_esati= zeros(sat_count, 3);
% position of all satellites in ECEF
est_vel_e_esati= zeros(sat_count, 3);
% Velocity of all satellites in ECEF
for index = 1: sat_count
    [est_pos_e_esati(index, :) , est_vel_e_esati(index, :)]= ...
        Satellite_position_and_velocity(time, sat_index(index));
end
while 1>0
% infinite loop unless be break out 
    pred_range_asati= zeros(sat_count, 1);
    % predicted pseudo ranges from current estimate of antenna position and
    % satellites 
    u_e_asati= zeros(sat_count, 3);
    % line-of-sight unit vector from current estimate of antenna position and
    % satellites
    for index = 1: numel(sat_index)
        [pred_range_asati(index, 1), u_e_asati(index, :)]= ...
            line_of_sight_vector(est_pos_e_esati(index, :)', pred_pos_e_ea);
    end
    pred_state= [pred_pos_e_ea; pred_clock_offset];
    % predicted state vector [ pos_x, pos_y pos_z clock offset]
    pred_meas_innov= obs_prange_a_sati- pred_range_asati - pred_clock_offset;
    % predicted measurement innovation vector
    meas_mat= horzcat(-u_e_asati, ones( sat_count, 1));
    % measurement matrix
    A= meas_mat'* meas_mat;
    B= inv(A)* meas_mat';
    C= B* pred_meas_innov;
    est_state= pred_state + C;
    % estimated state using unweighted least-squares
    x_change= abs(est_state(1)- pred_state(1));
    y_change= abs(est_state(2)- pred_state(2));
    z_change= abs(est_state(3)- pred_state(3));
    % absolute value of change in position of antenna in ECEF
    resultx = x_change < 0.01;
    resulty = y_change < 0.01;
    resultz = z_change < 0.01;
    % checking if absolute change in position in all axes is less than 10cm
    if resultx==1 && resulty==1 && resultz==1
    % If the estimated position solution has not moved by more than 10cm
    % in any axis we break out and return the current state estimates
        break
    end
    pred_pos_e_ea= est_state(1:3, 1);
    pred_clock_offset= est_state(4, 1);
    % setting estimated values as predicted values that will be used in
    % next iteration
end




