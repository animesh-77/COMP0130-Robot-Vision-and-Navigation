function [x_cap_new] = Single_Epoch_position(time_index, time, ...
    r_cap_e_ea_old, ~, filename)
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
% filename : name of the pseudo_range.csv file
% 
% Returns
% x_cap_new : position of the antenna




Define_Constants

delta_rho_a_c_old= 1;

data = readmatrix(filename);


%    Where i = index for satellite
% r_cap_e_esati_old = position for satellite i
%    Where i = index for satellite
% v_cap_e_esati_old = velocity for satellite i

sat_index= data(1, 2:end);

r_cap_e_esati_old= zeros(numel(sat_index), 3);
v_cap_e_esati_old= zeros(numel(sat_index), 3);

for index = 1: numel(sat_index)
    [r_cap_e_esati_old(index, :) , v_cap_e_esati_old(index, :)]= ...
        Satellite_position_and_velocity(time, sat_index(index));
end

% r_cap_e_esati_old is a 8x3 matrix


% rho_obs_sati_a where i is index of satellite
    
% Read the data using csvread

rho_obs_sati_a= data(time_index+1, 2:end)';
% only read the row at time step required


i= 0;
while 1>0
    
    
    % disp(r_cap_e_ea_old)
    %    Where i is the index for satellite
    % r_cap_e_asati_old = r_cap_e_esati_old -  r_cap_e_ea_old
    
    r_cap_e_asati_old= r_cap_e_esati_old- r_cap_e_ea_old';
    
    % magnitude of vector from antenna to satellite
    r_cap_e_asati_old_mag= sqrt(sum(power(r_cap_e_asati_old, 2), 2));
    
    %  where i is the index of the satellite
    % r_cap_e_asati_corr_old
    % 
    %  where i is the index of the satellite
    % u_e_asati
    
    r_cap_e_asati_corr_old= zeros(size(r_cap_e_asati_old_mag));
    
    u_e_asati= zeros(numel(r_cap_e_asati_old_mag), 3);
    

    C_I_e= [0 omega_ie/c 0; -omega_ie/c 0 0; 0 0 0];
    
    for index = 1: numel(sat_index)

        [r_cap_e_asati_corr_old(index, 1), u_e_asati(index, :)]= ...
            line_of_sight_vector(r_cap_e_esati_old(index, :)', r_cap_e_ea_old);

    end
    
    x_cap_old= [r_cap_e_ea_old; delta_rho_a_c_old];
    
    % fprintf('iteration %d x %-0.3f %-0.3f %-0.3f %-0.3f\n', i, x_cap_old);

    
    
    % delta_z_old
    
    delta_z_old= rho_obs_sati_a- r_cap_e_asati_corr_old - delta_rho_a_c_old;
    
    
    % H_e_G
    % 
    H_e_G= horzcat(-u_e_asati, ones(numel(sat_index), 1));
    
    A= H_e_G'* H_e_G;
    
    B= inv(A)* H_e_G';
    
    C= B* delta_z_old;
    
    x_cap_new= x_cap_old + C;
    
    result1 = (x_cap_new(1, 1)- x_cap_old(1, 1)) > .010;
    result2 = (x_cap_new(2, 1)- x_cap_old(2, 1)) > .010;
    result3 = (x_cap_new(3, 1)- x_cap_old(3, 1)) > .010;

    % if i < printi
    % 
    %     fprintf("%d %d %d\n", result1, result2, result3)
    % end


    if result1==0 && result2==0 && result3==0
        break
    end




    r_cap_e_ea_old= x_cap_new(1:3, 1);
    delta_rho_a_c_old= x_cap_new(4, 1);

    i= i+1;
    


end




