function updatedMeasurements = detectOutliers(measurements, delta_z, h)
    number_measurements = length(measurements);
    identity_matrix = eye(number_measurements);

    % Compute the residuals vector
    v = (h * inv(h' * h) * h' - identity_matrix) * delta_z;

    est = 5; % measurement error standard deviation
    est_squared = est^2;

    % Compute the Residuals Covariance Matrix
    c = (identity_matrix - h * inv(h' * h) * h') * est_squared;

    outlier_threshold = 6; % suggested value in Lab 1
    max_sat_range = max(measurements);
    min_sat_range = min(measurements);

    outlier_detected = false; 
    outlier_measurement_index = 0;
    max_residual = -inf;

    % Assuming c is at least an 8x8 matrix
    for i = 1:number_measurements
        % Normalize the residual
        norm_residual = (measurements(i) - min_sat_range) / (max_sat_range - min_sat_range);
        threshold = sqrt(c(i, i)) * outlier_threshold;
        
        % Check if there is an outlier
        if norm_residual > threshold
            fprintf('Measurement %f is an outlier with normalized residual %f\n', measurements(i), norm_residual);
            outlier_detected = true; % Set the flag
            if norm_residual > max_residual
                max_residual = norm_residual;
                outlier_measurement_index = i; % Store the index of the outlier
            end
        end
    end

    if outlier_detected
        % Remove the outlier measurement
        measurements(outlier_measurement_index) = []; 
    end

    % Return the updated measurements
    updatedMeasurements = measurements;
end
