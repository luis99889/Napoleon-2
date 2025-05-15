function [closest_sat_coords, closest_sat_indices, closest_sat_elevations_discrete, closest_sat_dists, num_times] = WP2_function(P, S, semiMajorAxis, inclination, gs_lat, gs_long, gs_alt, raan)

% Initializing the satellite scenario
startTime = datetime(2025,3,02,0,0,0);
stopTime = startTime + hours(2);
sampleTime = 30;
sc = satelliteScenario(startTime,stopTime,sampleTime);

% Create an array for time samples in date format
times = startTime : seconds(sampleTime) : stopTime;

% Total number of satellites    
total_sats = P * S;

% Define constant orbital parameters for all satellites
eccentricity = zeros(1, total_sats);           % Circular orbits
argumentOfPerigee = zeros(1, total_sats);      % Arbitrary for circular orbits
semiMajorAxis = repmat(semiMajorAxis, 1, total_sats); % Same for all
inclination = repmat(inclination, 1, total_sats);     % Same for all

% Compute RAAN: repeat each plane's RAAN for S satellites
RAAN = repelem((0:P-1) * raan / P, S);

% Compute true anomaly: repeat the sequence [0, 360/S, ..., (S-1)*360/S] for P planes
trueAnomaly = repmat((0:S-1) * 360 / S, 1, P);

% Create the satellite array
sat = satellite(sc, semiMajorAxis, eccentricity, inclination, RAAN, ...
                    argumentOfPerigee, trueAnomaly);

% Define the GS (GU) in the Geodetic Coordinates
gs = groundStation(sc, gs_lat, gs_long);

% Compute GS ECEF position (no reference model specified, defaults to spherical Earth)
gs_ecef_pos = lla2ecef([gs_lat, gs_long, gs_alt]);
gs_ecef_pos_time = repmat(gs_ecef_pos', 1, numel(times));

% Define dimensions
num_times = numel(times);
num_sats = numel(sat);

% Define array to store positions of satellites through time (3D: 3 coords × times × satellites)
all_sat_pos = states(sat, "CoordinateFrame", "ecef");

% Define array for Access status to signal visibility between GS and sats
ac = access(sat, gs);
intvls = accessIntervals(ac);

% Initialize status array to store logical access status at each time
status = false(num_times, num_sats);  % Preallocate with logical false

% Extract status for each time step
for j = 1:num_times
    status(j, :) = accessStatus(ac, times(j));  % Store logical status for each satellite
end

% Define the number of closest satellites to find
n = 3;

% Preallocate matrices to store data for n closest satellites
closest_sat_coords = zeros(num_times, n, 3);  % 3D matrix: num_times x n x 3 (for x,y,z)
closest_sat_indices = zeros(num_times, n);    % 2D matrix: num_times x n (for satellite indices)
closest_sat_dists = zeros(num_times, n);      % 2D matrix: num_times x n (for distances in meters)
closest_sat_elevations = zeros(num_times, n); % 2D matrix: num_times x n (for elevation angles in degrees)
closest_sat_elevations_discrete = zeros(num_times, n); % 2D matrix: num_times x n (for discretized elevation angles)

% First Loop: Find the n closest satellites at each time and store their data (coords, indices, distances)
fprintf('Calculating closest satellites and distances...\n'); % Progress indicator
for j = 1:numel(times)
    active_sats = find(status(j, :));
    gs_pos = gs_ecef_pos_time(:, j); % GS ECEF position at current time

    if ~isempty(active_sats)
        % Calculate distances only for active satellites
        distances = zeros(1, length(active_sats));
        sat_positions_active = zeros(3, length(active_sats)); % Store positions temporarily

        for k = 1:length(active_sats)
            sat_idx = active_sats(k);
            sat_pos = all_sat_pos(:, j, sat_idx);
            sat_positions_active(:, k) = sat_pos; % Store position for later use
            distances(k) = norm(gs_pos - sat_pos);
        end

        % Sort distances to find the n closest satellites
        [sorted_dists, dist_idx] = sort(distances);

        % Limit to n or the number of active satellites, whichever is smaller
        num_closest = min(n, length(active_sats));
        closest_sats_indices_in_active = dist_idx(1:num_closest); % Indices within the 'active_sats' array
        closest_sats = active_sats(closest_sats_indices_in_active); % Actual satellite indices (1 to num_sats)
        closest_dists = sorted_dists(1:num_closest);

        % Store ECEF coordinates, indices, and distances
        for k = 1:num_closest
            sat_idx = closest_sats(k); % The actual satellite index
            sat_pos = sat_positions_active(:, closest_sats_indices_in_active(k)); % Retrieve stored position

            % Store existing data
            closest_sat_coords(j, k, :) = sat_pos;          % Coordinates
            closest_sat_indices(j, k) = sat_idx;            % Satellite index
            closest_sat_dists(j, k) = closest_dists(k);     % Distance in meters
        end
    else
        % No satellites in view, matrices remain zeros for this time step
    end

    % Optional: Add progress update to console
    if mod(j, round(num_times/10)) == 0 || j == num_times
        fprintf('Processed %d/%d time steps (%.0f%%)\n', j, num_times, (j/num_times)*100);
    end
end

% Second Loop: Compute and store elevation angles for the n closest satellites using the geometrical model
fprintf('Calculating elevation angles for closest satellites...\n');
R_E = 6371000; % Earth's mean radius in meters (spherical assumption)
h_0 = semiMajorAxis(1) - R_E; % Satellite altitude in meters (fixed for circular orbit)
gs_lat_rad = deg2rad(gs_lat); % GS latitude in radians
gs_long_rad = deg2rad(gs_long); % GS longitude in radians

for j = 1:numel(times)
    % Get the indices of the n closest satellites for this time step
    closest_sats = closest_sat_indices(j, :);
    num_closest = sum(closest_sats > 0); % Number of valid satellites (non-zero indices)

    if num_closest > 0
        for k = 1:num_closest
            sat_idx = closest_sats(k);
            if sat_idx > 0 % Check if valid satellite index
                % Compute SSP coordinates using ecef2lla
                sat_pos = squeeze(closest_sat_coords(j, k, :)); % ECEF position of the satellite
                ssp_lla = ecef2lla(sat_pos'); % Convert ECEF to LLA
                ssp_lla(3) = 0; % Set SSP altitude to zero (sub-satellite point on Earth's surface)

                % Extract SSP latitude and longitude in radians
                ssp_lat_rad = deg2rad(ssp_lla(1));
                ssp_long_rad = deg2rad(ssp_lla(2));

                % Compute central angle gamma_ell using spherical law of cosines
                delta_phi = gs_long_rad - ssp_long_rad; % Longitude difference
                gamma_ell = acos(sin(ssp_lat_rad) * sin(gs_lat_rad) + ...
                                 cos(ssp_lat_rad) * cos(gs_lat_rad) * cos(delta_phi));

                % Compute elevation angle delta_ell using the formula from the paper
                if sin(gamma_ell) > eps % Avoid division by zero
                    term = (cos(gamma_ell) - (R_E / (R_E + h_0))) / sin(gamma_ell);
                    elevation = rad2deg(atan(term)); % Convert to degrees
                else
                    elevation = 90; % Satellite directly overhead
                end

                % Store elevation
                closest_sat_elevations(j, k) = elevation;
            end
        end
    end

    % Optional: Add progress update to console for elevation calculation
    if mod(j, round(num_times/10)) == 0 || j == num_times
        fprintf('Elevation processed %d/%d time steps (%.0f%%)\n', j, num_times, (j/num_times)*100);
    end
end

% Define the minimum elevation threshold
min_elevation = 15; % Example: degrees

% Discretize elevation angles into predefined bins [20, 30, 45, 60, 70]
discretized_values = [20 30 45 60 70];

% Replace the current discretization loop
for j = 1:num_times
    for k = 1:n
        elevation = closest_sat_elevations(j, k);
        
        if elevation > min_elevation % Only discretize angles above the threshold
            % Find the closest value in discretized_values
            [~, closest_idx] = min(abs(discretized_values - elevation));
            closest_sat_elevations_discrete(j, k) = discretized_values(closest_idx);
        else
            % Set angles below the threshold to 0
            closest_sat_elevations_discrete(j, k) = 0;
        end
    end
end

fprintf('Discretization complete.\n');

% Visualize (optional)
play(sc);

% Highlight ground station
gs.MarkerColor = [1 0 0]; % Red marker for ground station

end
