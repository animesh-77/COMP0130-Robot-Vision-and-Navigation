function [est_state]= Single_epoch_velocity(time_index, time, ...
    pred_r_ea, range_rate_filename)


Define_Constants;

% Read the psudo_range_rate.csv file
data = readmatrix(range_rate_filename);

meas_range_rate_asati= data(time_index+1, 2:end)'; %(satellite count, 1)

% Getting position and velocity of all satellites

sat_index= data(1, 2:end);
sat_count= numel(sat_index);

pred_r_esati= zeros(sat_count, 3);
pred_v_esati= zeros(sat_count, 3);

for index = 1: sat_count
    [pred_r_esati(index, :) , pred_v_esati(index, :)]= ...
    Satellite_position_and_velocity(time, sat_index(index));
end


% starting with zero velocities
pred_v_ea_resolved= zeros(3, 1);

pred_v_ea=  zeros(1, 3);  
clock_drift= 0;

pred_state= zeros(4, 1);
pred_state(1:3, 1)= pred_v_ea';
pred_state(4, 1)= clock_drift;

pred_range_rate_asati= zeros(8, 1);

i= 0;

while i>=0
    
    % fprintf("iter %d - ", i)
    % fprintf("N : %.2f m/s E : %.2f m/s D : %.2f m/s\n", pred_v_ea_resolved)

    pred_range_asati= zeros(sat_count, 1);
    pred_u_asati= zeros(sat_count, 3);
    
    for index=1: sat_count

       [pred_range_asati(index, 1), pred_u_asati(index, :)]= ...
           line_of_sight_vector(pred_r_esati(index, :)', pred_r_ea);
    end
    




    for index = 1: sat_count
    
        C_I_e= eye(3); % (3, 3)
        C_I_e(1, 2)= omega_ie*pred_range_asati(index, 1)/c;
        C_I_e(2, 1)= -omega_ie*pred_range_asati(index, 1)/c;
        
        A= pred_v_esati(index, :)'+ Omega_ie*pred_r_esati(index, :)'; % (3, 1)
    
        B= pred_v_ea'+ Omega_ie* pred_r_ea; % (3, 1)
    
    
        pred_range_rate_asati(index, 1)= pred_u_asati(index, :)*((C_I_e* A)- B);
    

    end

    pred_meaus_innov= meas_range_rate_asati- pred_range_rate_asati- pred_state(4, 1);

    H= ones(8, 4);
    H(:, 1:3)= -pred_u_asati;

    est_state= pred_state+ inv(H'*H)*H'*pred_meaus_innov;


    [~, ~, ~, est_v_ea_resolved]= pv_ECEF_to_NED(pred_r_ea, est_state(1:3, 1));

    Delta_vx = est_v_ea_resolved(1, 1)- pred_v_ea_resolved(1, 1);
    Delta_vy = est_v_ea_resolved(2, 1)- pred_v_ea_resolved(2, 1);
    Delta_vz = est_v_ea_resolved(3, ...
        1)- pred_v_ea_resolved(3, 1);
    
    if Delta_vx < 0.01 && Delta_vy < 0.01 && Delta_vz < 0.01
        break
    end


    i = i+1;

    pred_v_ea=  est_state(1:3, 1)'; 
    clock_drift= est_state(4, 1);
    pred_state= est_state;
    pred_v_ea_resolved= est_v_ea_resolved;

    % break

end


end