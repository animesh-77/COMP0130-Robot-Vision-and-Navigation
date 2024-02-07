function [est_state] = Single_Epoch_position(time, ...
    pred_r_e_ea, ~, sat_index, obs_prange_a_sati)
% 
% 
% Function to get single epoch position solution without any priors
% Iterates until position solution does not improve by more than 10cm in any axis
% 
% Inputs :
% time_index : index of time in the pseudo_range.csv file (1, 1)
% time : time in seconds (1, 1)
% r_cap_e_ea_old : prior position of the antenna. Can also be centre of earth (3, 1)  
% ~ : velocity vector not required . (3, 1)
% sat_index : index number of satellites 
% obs_prange_a_sati : observed pseudo range measurements from all satellites
% 
% Returns
% x_cap_new : position of the antenna




Define_Constants

pred_clock_offset= 1;

% data = readmatrix(filename);


%    Where i = index for satellite
% r_cap_e_esati_old = position for satellite i
%    Where i = index for satellite
% v_cap_e_esati_old = velocity for satellite i

sat_count= numel(sat_index);

est_r_e_esati= zeros(sat_count, 3);
est_v_e_esati= zeros(sat_count, 3);

for index = 1: sat_count
    [est_r_e_esati(index, :) , est_v_e_esati(index, :)]= ...
        Satellite_position_and_velocity(time, sat_index(index));
end

% est_r_e_esati is a number of satellites x 3 matrix


% rho_obs_sati_a where i is index of satellite
    
% obs_prange_a_sati = data(time_index+1, 2:end)';
% only read the row at time step required


i= 0;
while 1>0
        
    %  where i is the index of the satellite
    % r_cap_e_asati_corr_old
    % 
    %  where i is the index of the satellite
    % u_e_asati
    
    pred_r_e_asati_corr= zeros(sat_count, 1);
    
    u_e_asati= zeros(sat_count, 3);
    
    
    for index = 1: numel(sat_index)

        [pred_r_e_asati_corr(index, 1), u_e_asati(index, :)]= ...
            line_of_sight_vector(est_r_e_esati(index, :)', pred_r_e_ea);

    end
    
    pred_state= [pred_r_e_ea; pred_clock_offset];
    
    % fprintf('iteration %d x %-0.3f %-0.3f %-0.3f %-0.3f\n', i, x_cap_old);

    
    
    % delta_z_old
    
    pred_meas_innov= obs_prange_a_sati- pred_r_e_asati_corr - pred_clock_offset;
    
    
    % H_e_G
    % 
    H_e_G= horzcat(-u_e_asati, ones( sat_count, 1));
    
    A= H_e_G'* H_e_G;
    
    B= inv(A)* H_e_G';
    
    C= B* pred_meas_innov;
    
    est_state= pred_state + C;
    

    x_change= abs(est_state(1)- pred_state(1));
    y_change= abs(est_state(2)- pred_state(2));
    z_change= abs(est_state(3)- pred_state(3));

    resultx = x_change < 0.01;
    resulty = y_change < 0.01;
    resultz = z_change < 0.01;

    % if i < printi
    % 
    %     fprintf("%d %d %d\n", result1, result2, result3)
    % end


    if resultx==1 && resulty==1 && resultz==1
    % if result4 == 1
        % [m, n] = size(est_state);
        break
    end




    pred_r_e_ea= est_state(1:3, 1);
    pred_clock_offset= est_state(4, 1);

    i= i+1;
    


end




