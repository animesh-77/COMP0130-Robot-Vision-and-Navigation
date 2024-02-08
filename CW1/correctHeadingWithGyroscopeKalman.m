function correctedHeading = correctHeadingWithGyroscopeKalman(data)
    % Extract relevant data
    gyroRates = data(:, 6);  % Gyroscope angular rate in rad/s
    compassHeadings = deg2rad(data(:, 7)); % Compass headings in radians
    
    % Constants
    deltaTime = 0.5;  % Time step in seconds, adjust based on your data sampling rate
    gyro_error_variance = 3 * 10^(-6);  % Gyro measurement error variance in rad^2/s^3
    compass_noise_variance = 10^(-8); % Compass measurement noise variance in rad^2
    
    % Initialize state and covariance matrix
    states = [compassHeadings(1); 0]; % Initial heading and rate of change of heading (rad/s)
    P = [compass_noise_variance, 0; 
         0, gyro_error_variance]; % Initial error covariance matrix
    
    % System noise covariance matrix (Q)
    Q = [1/4 * deltaTime^4, 1/2 * deltaTime^3; 
         1/2 * deltaTime^3, deltaTime^2] * gyro_error_variance;
    
    % Measurement matrix (H) for heading update
    H = [1, 0];
    
    % Measurement noise covariance matrix (R)
    R = compass_noise_variance;
    
    % Allocate space for correctedHeading
    correctedHeading = zeros(size(compassHeadings));
    correctedHeading(1) = compassHeadings(1); % Set initial corrected heading
    
    for i = 2:length(compassHeadings)
        % Prediction using gyro rates
        phi = [1, deltaTime; 0, 1]; % State transition model
        gyroRate = gyroRates(i-1); % Gyro rate at previous timestep
        predictedRateOfChange = gyroRate; % Assuming gyroRate is the rate of change of heading
        states = phi * states + [0; predictedRateOfChange]; % Predict next state with gyro adjustment
        P = phi * P * phi' + Q; % Predict next covariance
        
        % Update using compass heading
        K = P * H' / (H * P * H' + R); % Kalman gain
        z = compassHeadings(i); % New measurement from compass
        y = z - H * states; % Measurement residual
        states = states + K * y; % Update state estimate with compass measurement
        P = (eye(2) - K * H) * P; % Update covariance estimate
        
        correctedHeading(i) = states(1); % Store corrected heading
    end
    
    % Convert corrected headings back to degrees
    correctedHeading = rad2deg(correctedHeading);
end
