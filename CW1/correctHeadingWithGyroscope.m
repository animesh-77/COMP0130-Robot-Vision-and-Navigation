function correctedHeading = correctHeadingWithGyroscope(data)
    % Extract gyroscope data and compass headings
    gyroRates = data(:, 6); % Gyroscope angular rates in rad/s
    compassHeadings = deg2rad(data(:, 7)); % Convert compass headings to radians
    
    deltaTime = 0.5;
    correctedHeading = zeros(size(compassHeadings));
    
    % Assume initial heading is correct
    correctedHeading(1) = compassHeadings(1);
    
    % Integrate gyroscope data to correct heading
    for i = 2:length(compassHeadings)
        % Integrate gyroscope data to estimate change in heading
        deltaHeading = gyroRates(i-1) * deltaTime;
        
        % Correct the heading
        correctedHeading(i) = correctedHeading(i-1) + deltaHeading;
        
        % Fuse with compass heading with simple averaging
        correctedHeading(i) = mean([correctedHeading(i), compassHeadings(i)]);
    end
    
    % Convert back to degrees
    correctedHeading = rad2deg(correctedHeading);
end
