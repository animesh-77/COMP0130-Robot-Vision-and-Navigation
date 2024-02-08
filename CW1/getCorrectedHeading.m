data = readmatrix("Dead_reckoning.csv");

correctedHeading = correctHeadingWithGyroscope(data);

% originalHeading = data(:, 7); % Extract original heading in degrees

% Combine original and corrected headings into a two-column matrix
% combinedHeadings = [originalHeading, correctedHeading];

% Display the combined headings
% disp('Original Heading (deg) | Corrected Heading (deg)');
% disp(combinedHeadings);

data(:, 7) = correctedHeading;

% disp(data)
