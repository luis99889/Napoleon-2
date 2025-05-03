clc;
clear;
close all;

% Initializing the satellite scenario
startTime = datetime(2025,3,02,0,0,0);
stopTime = startTime + hours(2);
sampleTime = 30;
sc = satelliteScenario(startTime,stopTime,sampleTime);

% Create an array for time samples in date format
times = startTime : seconds(sampleTime) : stopTime;

% Define constellation parameters
P = 6;                % Number of orbit planes
S = 11;               % Satellites per orbit
semiMajorAxis = 6371000 + 780000; % Orbital radius in meters (R_E + altitude)
inclination = 86.4;   % degrees

% Add satellites using the constellation function
sat = createConstellation(sc, P, S, semiMajorAxis, inclination);

% Set Turin's coordinates for the fixed GU
gs_lat = 45.07;
gs_long = 7.69;
gs_alt = 0;

% Define the GS (GU) in the Geodetic Coordinates
gs = groundStation(sc, gs_lat, gs_long);

% Compute GS ECEF position (no reference model specified, defaults to spherical Earth)
gs_ecef_pos = lla2ecef([gs_lat, gs_long, gs_alt]);
gs_ecef_pos_time = repmat(gs_ecef_pos', 1, numel(times));

% % Verify computed GS radius to confirm spherical Earth assumption
% R_E_computed = norm(gs_ecef_pos);
% fprintf('Computed GS radius: %.2f meters\n', R_E_computed);
% 
% % Debug: Print GS coordinates for verification
% fprintf('Ground Station Coordinates: Lat = %.4f deg, Lon = %.4f deg, Alt = %.4f m\n', gs_lat, gs_long, gs_alt);
% fprintf('Ground Station ECEF Position: [%.2f, %.2f, %.2f] m\n', gs_ecef_pos(1), gs_ecef_pos(2), gs_ecef_pos(3));
% 
% % Convert back to LLA to confirm consistency
% LLA = ecef2lla([gs_ecef_pos(1) gs_ecef_pos(2) gs_ecef_pos(3)]);
% fprintf('Ground Station ECEF to LLA Check: Lat = %.4f deg, Lon = %.4f deg, Alt = %.4f m\n', LLA(1), LLA(2), LLA(3));

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
h_0 = semiMajorAxis - R_E; % Satellite altitude in meters (fixed for circular orbit)
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

% % Third Loop: Discretize elevation values into 5 bins (0 to 90 degrees)
% fprintf('Discretizing elevation angles into 5 bins...\n');
% bin_size = 90 / 5; % 18 degrees per bin
% bin_boundaries = 0:bin_size:90; % Boundaries: [0, 18, 36, 54, 72, 90]
% bin_centers = (bin_boundaries(1:end-1) + bin_boundaries(2:end)) / 2; % Centers: [9, 27, 45, 63, 81]
% 
% for j = 1:num_times
%     for k = 1:n
%         elevation = closest_sat_elevations(j, k);
%         if elevation > 0 % Only discretize non-zero (valid) elevations ..... 20 30 45 60 70 function quantize
%             % Find the appropriate bin using boundaries
%             if elevation < bin_boundaries(2)
%                 closest_sat_elevations_discrete(j, k) = bin_centers(1); % Bin 1: 9 degrees
%             elseif elevation < bin_boundaries(3)
%                 closest_sat_elevations_discrete(j, k) = bin_centers(2); % Bin 2: 27 degrees
%             elseif elevation < bin_boundaries(4)
%                 closest_sat_elevations_discrete(j, k) = bin_centers(3); % Bin 3: 45 degrees
%             elseif elevation < bin_boundaries(5)
%                 closest_sat_elevations_discrete(j, k) = bin_centers(4); % Bin 4: 63 degrees
%             else
%                 closest_sat_elevations_discrete(j, k) = bin_centers(5); % Bin 5: 81 degrees
%             end
%         else
%             closest_sat_elevations_discrete(j, k) = 0; % Keep invalid entries as 0
%         end
%     end
% end

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

% Optional: Display a sample of the stored data for verification
time_step_to_display = 60; % Example time step
fprintf('\nSample data for n=%d closest satellites at time step %d:\n', n, time_step_to_display);
fprintf('Satellite Indices:\n');
disp(closest_sat_indices(time_step_to_display, :));
fprintf('Distances (km):\n');
disp(closest_sat_dists(time_step_to_display, :) / 1000); % Convert to km
fprintf('Original Elevation Angles (degrees):\n');
disp(closest_sat_elevations(time_step_to_display, :));
fprintf('Discretized Elevation Angles (degrees):\n');
disp(closest_sat_elevations_discrete(time_step_to_display, :));
fprintf('ECEF Coordinates (m):\n');
disp(squeeze(closest_sat_coords(time_step_to_display, :, :))); % Display as n x 3 matrix

% Visualize (optional)
play(sc);

% Highlight ground station
gs.MarkerColor = [1 0 0]; % Red marker for ground station

% Function to create constellation
function sat = createConstellation(sc, P, S, semiMajorAxis, inclination)
    % Total number of satellites
    total_sats = P * S;

    % Define constant orbital parameters for all satellites
    eccentricity = zeros(1, total_sats);           % Circular orbits
    argumentOfPerigee = zeros(1, total_sats);      % Arbitrary for circular orbits
    semiMajorAxis = repmat(semiMajorAxis, 1, total_sats); % Same for all
    inclination = repmat(inclination, 1, total_sats);     % Same for all

    % Compute RAAN: repeat each plane's RAAN for S satellites
    RAAN = repelem((0:P-1) * 180 / P, S);

    % Compute true anomaly: repeat the sequence [0, 360/S, ..., (S-1)*360/S] for P planes
    trueAnomaly = repmat((0:S-1) * 360 / S, 1, P);

    % Create the satellite array
    sat = satellite(sc, semiMajorAxis, eccentricity, inclination, RAAN, ...
                    argumentOfPerigee, trueAnomaly);
end
%% WP3

% Parameters
seed = 72;
Sample_Rate = 400; % Hz
numSamples = Sample_Rate * num_times * 60; % Sample rate times duration of simulation in seconds

all_faded_waves = cell(num_times, 1);
all_channel_gains = cell(num_times, 1);
all_states = cell(num_times, 1);
all_distances = cell(num_times, 1);

Sat_Ang_time = zeros(num_times,1);

% Loop over time steps
new_State = "good";  % The very first state is set 'good' (lowercase)

for t = 1:num_times
    % Use elevation for this time
    elev = closest_sat_elevations_discrete(t);
    
    Sat_Ang_time(t) = elev;
    
    if elev > 0  % Only simulate if elevation is valid
        chan = p681LMSChannel;
        chan.Environment = "Urban";
        chan.CarrierFrequency = 3.8e9;
        chan.ElevationAngle = elev;
        chan.MobileSpeed = 0.5;
        chan.AzimuthOrientation = 0;
        chan.SampleRate = Sample_Rate;
        
        chan.InitialState = new_State;  % Use the last state of the previous iteration
        
        chan.FadingTechnique = "Filtered Gaussian noise";
        chan.RandomStream = "mt19937ar with seed";
        chan.Seed = seed + t; % Ensure unique seed for each time

        % Set random number generator with seed
        rng(seed);
        
        % Channel duration 60 sec, because we update the angle every 60 sec
        chanDur = 60; 
        % Random input waveform
        numSamples = floor(chan.SampleRate * chanDur) + 1;
        in = complex(randn(numSamples,1), randn(numSamples,1));
        
        % Pass the input signal through channel
        [fadedWave, channelCoefficients, sampleTimes, stateSeries] = step(chan, in);

        % Store results
        all_faded_waves{t} = fadedWave;
        all_channel_gains{t} = channelCoefficients;
        all_states{t} = stateSeries;
        all_times{t} = sampleTimes;

        % Memorize the last state, converted to lowercase
        previousState = stateSeries(end);

        if previousState == 0

            new_State = "bad";

        else 
            new_State = "good";
        end

                  
    else
        % No satellite in view — store NaN
        all_faded_waves{t} = NaN;
        all_channel_gains{t} = NaN;
        all_states{t} = NaN;
        all_times{t} = NaN;        
    end
end



% Total number of simulations
num_total = numel(all_faded_waves);

% Number of simulations to concatenate
num_to_concat = 20;

% Indices of the simulations to be concatenated
selected_idxs = num_total - num_to_concat + 1 : num_total;

% Initialize arrays for the concatenated data
full_time = [];
full_fadedWave = [];
full_inputWave = [];
full_channel = [];
full_state = [];

t_offset = 0; % Cumulative time offset

% Loop over the selected simulations
for idx = selected_idxs
    fadedWave = all_faded_waves{idx};
    channelCoefficients = all_channel_gains{idx};
    stateSeries = all_states{idx};
    timeVector = all_times{idx};
    
    % Align the time vector with the cumulative time offset
    adjusted_time = timeVector + t_offset;
    
    % Append the current simulation data to the concatenated arrays
    full_time = [full_time; adjusted_time(:)];
    full_fadedWave = [full_fadedWave; fadedWave(:)];
    full_channel = [full_channel; channelCoefficients(:)];
    full_state = [full_state; stateSeries(:)];
    
    % Assume the input waveform "in" is the same for all segments
    inputSegment = in(1:length(timeVector));
    full_inputWave = [full_inputWave; inputSegment(:)];
    
    % Update the cumulative time offset
    t_offset = t_offset + timeVector(end);
end


Last_20_Ang = Sat_Ang_time(end-19:end);


fprintf('Last 20 Elevation Angles: ');
fprintf('%f ', Last_20_Ang);
fprintf('\n');

% Power Profile
figure(1)
plot(full_time, 20*log10(abs(full_inputWave)), ...
     full_time, 20*log10(abs(full_fadedWave)))
title('Power Profile - Last 20 Simulations')
legend('Input Waveform', 'Faded Waveform')
xlabel('Time (s)')
ylabel('Power (dB)')

% Channel Gain
figure(2)
plot(full_time, 20*log10(abs(full_channel)))
title('Channel Gain - Last 20 Simulations')
xlabel('Time (s)')
ylabel('Path Gain (dB)')

% State Series 
figure(3)
plot(full_time, full_state)
title('State Series - Last 20 Simulations')
axis([0 full_time(end) -0.5 1.5])
xlabel('Time (s)')
ylabel('State')

%%
figure(4)
plot(Sat_Ang_time, 'o-', 'MarkerSize', 3, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
title('Satellite Elevation Angle over Time')
xlabel('Time Step')
ylabel('Elevation Angle (degrees)')

figure(5)
histogram(Sat_Ang_time(~isnan(Sat_Ang_time)), 'BinWidth', 5, 'FaceColor', 'g')
title('Elevation angles occurrency')
xlabel('Elevazion (degrees)')
ylabel('Occurrencies')
grid on

%% WP4

% free space path loss (distance), tx power, gain of antenna, bandwidth, noise psd

% maximal ratio combining