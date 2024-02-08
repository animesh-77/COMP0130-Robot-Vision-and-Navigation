function [est_state]= Single_epoch_velocity(time, ...
    est_pos_e_ea, sat_index, osb_range_rate_asati)
% 
% Function to estimate velocity of antenna using unweighted
% least squares estimation
% 
% Inputs 
% time : time in seconds (1, 1)
% est_pos_e_ea : estimated position on antenna in ECEF (3, 1)
% sat_index :index number of satellites (n,) where n is the number of satellites
% osb_range_rate_asati : observed pseudo range rate measurements from all satellites (n, 1)
% 
% Returns
% est_state : velocity of the antenna and clock drift
% 
Define_Constants;
% Getting position and velocity of all satellites
sat_count= numel(sat_index);
pred_r_esati= zeros(sat_count, 3);
% position of all satellites in ECEF
pred_v_esati= zeros(sat_count, 3);
% velocity of all satellites in ECEF
for index = 1: sat_count
    [pred_r_esati(index, :) , pred_v_esati(index, :)]= ...
    Satellite_position_and_velocity(time, sat_index(index));
end
% starting with zero velocities for the antenna
pred_v_ea=  zeros(1, 3);  
clock_drift= 0;
% intial velocity in NED frame
[~, ~, ~, pred_vel_ea_resolved]= pv_ECEF_to_NED(est_pos_e_ea, pred_v_ea');
% predicted state vector
pred_state= zeros(4, 1);
pred_state(1:3, 1)= pred_v_ea';
pred_state(4, 1)= clock_drift;
pred_range_rate_asati= zeros(8, 1);
while 1>0

    pred_range_asati= zeros(sat_count, 1);
    % predicted pseudo ranges from current estimate of antenna position and
    % satellites 
    pred_u_asati= zeros(sat_count, 3);
    % line-of-sight unit vector from current estimate of antenna position and
    % satellites
    for index=1: sat_count
       [pred_range_asati(index, 1), pred_u_asati(index, :)]= ...
           line_of_sight_vector(pred_r_esati(index, :)', est_pos_e_ea);
    end
    for index = 1: sat_count
        % iterating over all satellites
        C_I_e= eye(3); % (3, 3)
        C_I_e(1, 2)= omega_ie*pred_range_asati(index, 1)/c;
        C_I_e(2, 1)= -omega_ie*pred_range_asati(index, 1)/c;
        A= pred_v_esati(index, :)'+ Omega_ie*pred_r_esati(index, :)'; % (3, 1)
        B= pred_v_ea'+ Omega_ie* est_pos_e_ea; % (3, 1)
        pred_range_rate_asati(index, 1)= pred_u_asati(index, :)*((C_I_e* A)- B);
        % predicted range rate measurement for a satellite
    end
    % predicted measurement innovation vector
    pred_meaus_innov= osb_range_rate_asati- pred_range_rate_asati- pred_state(4, 1);
    % measurement matrix
    H= ones(8, 4);
    H(:, 1:3)= -pred_u_asati;
    % estimated state using unweighted least-squares
    est_state= pred_state+ inv(H'*H)*H'*pred_meaus_innov;
    % resolving velocity in NED frame
    [~, ~, ~, est_v_ea_resolved]= pv_ECEF_to_NED(est_pos_e_ea, est_state(1:3, 1));
    % absolute value of change in velocity of antenna in NED
    Delta_vx = abs(est_v_ea_resolved(1, 1)- pred_vel_ea_resolved(1, 1));
    Delta_vy = abs(est_v_ea_resolved(2, 1)- pred_vel_ea_resolved(2, 1));
    Delta_vz = abs(est_v_ea_resolved(3, ...
        1)- pred_vel_ea_resolved(3, 1));
    % checking if absolute change in velocity in all axes is less than 10cm/s
    if Delta_vx < 0.01 && Delta_vy < 0.01 && Delta_vz < 0.01
        break
    end
    % setting estimated values as predicted values that will be used in
    % next iteration
    pred_v_ea=  est_state(1:3, 1)'; 
    clock_drift= est_state(4, 1);
    pred_state= est_state;
    pred_vel_ea_resolved= est_v_ea_resolved;
end
end