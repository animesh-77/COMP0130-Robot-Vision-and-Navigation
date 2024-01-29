function [pred_range_asati, pred_u_asati] = line_of_sight_vector(pred_r_esati, ...
    pred_r_ea)
% Calculate range and line of sight vector unit vector
% 
% 
% Input :
% pred_r_esati predicted r vector from center of earth to satellite (3, 1)
% pred_r_ea predicted r vector from centre of earth to antenna (3, 1)
% 
% 
% Returns : 
% pred_range_asati - predicted range from antenna to satellite (1, 1)
% pred_u_asati - line-of-sight unit vector (3, 1)
% 
% 


c = 299792458; % Speed of light in m/s
omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s


C= eye(3);

pred_range_asati = sqrt((C*pred_r_esati - pred_r_ea)' * (C*pred_r_esati - pred_r_ea)); % (1, 1)

C(2, 1)= -omega_ie* pred_range_asati/c;
C(1, 2)= omega_ie* pred_range_asati/c;

pred_range_asati = sqrt((C*pred_r_esati - pred_r_ea)' * (C*pred_r_esati - pred_r_ea)); % (1, 1)


pred_u_asati= (C*pred_r_esati - pred_r_ea)/pred_range_asati ; %(3, 1) matrix

end
