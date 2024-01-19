function [lat_deg, long_deg, height_m] = First_Epoch_pos_from_centre(time_index, time, ...
    r_cap_e_ea_old, v_cap_e_ea_old, ref_vec)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% get position and velocity from previous position


Define_Constants

delta_rho_a_c_old= 1;

%    Where i = index for satellite
% r_cap_e_esati_old = position for satellite i
%    Where i = index for satellite
% v_cap_e_esati_old = velocity for satellite i

sat_index= [2 17 18 22 23 26 27 28];

r_cap_e_esati_old= zeros(numel(sat_index), 3);
v_cap_e_esati_old= zeros(numel(sat_index), 3);

for index = 1: numel(sat_index)
    [r_cap_e_esati_old(index, :) , v_cap_e_esati_old(index, :)]= ...
        Satellite_position_and_velocity(time, sat_index(index));
end
% r_cap_e_esati_old is a 8x3 matrix


% rho_obs_sati_a where i is index of satellite
    
filename = 'Workshop1_Pseudo_ranges.csv';
% Read the data using csvread
data = readmatrix(filename);
rho_obs_sati_a= data(time_index+1, 2:end)';
% only read the row at time step required


i= 0;
while 1>0
    i= i+1;
    
    
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
    
        A= (C_I_e* r_cap_e_asati_old_mag(index, 1)+ diag(ones(1, 3)))* r_cap_e_esati_old(index, :)';
        
        B= A- r_cap_e_ea_old;
        
        C= sqrt(B'*B);
    
        r_cap_e_asati_corr_old(index, 1)= C;
    
        u_e_asati(index, :) = B'/C;
        
    end
    
    x_cap_old= [r_cap_e_ea_old; delta_rho_a_c_old];
    
    
    
    
    % delta_z_old
    
    delta_z_old= rho_obs_sati_a- r_cap_e_asati_corr_old - delta_rho_a_c_old;
    
    
    % H_e_G
    % 
    H_e_G= horzcat(-u_e_asati, ones(numel(sat_index), 1));
    
    A= inv(H_e_G'* H_e_G);
    
    B= A* H_e_G';
    
    C= B* delta_z_old;
    
    x_cap_new= x_cap_old + C;
    
    fprintf("%d", i)
    fprintf("%0.2f\n", x_cap_new(1:3, 1) )
    % if all( (x_cap_new(1:3, 1)- ref_vec ) < 10)
    %     break
    % end
    
    r_cap_e_ea_old= x_cap_new(1:3, 1);
    delta_rho_a_c_old= x_cap_new(4, 1);

    
    if i ==4
        break
    end


end

[lat_rad, long_rad, height_m]= pv_ECEF_to_NED(x_cap_new(1:3), [0 0 0]');

lat_deg= rad_to_deg*lat_rad;

long_deg= rad_to_deg* long_rad;


