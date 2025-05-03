function [closest_sat_coords, closest_sat_indices, closest_sat_elevations_discrete, closest_sat_dists, num_times] = WP2(P, S, semiMajorAxis, inclination, gs_lat, gs_long)

    % Initializing the satellite scenario
    startTime = datetime(2025,3,02,0,0,0);
    stopTime = startTime + hours(2);
    sampleTime = 30;
    sc = satelliteScenario(startTime,stopTime,sampleTime);

    % Create an array for time samples in date format
    times = startTime : seconds(sampleTime) : stopTime;
    num_times = numel(times); % Number of time steps

    % Total number of satellites
    total_sats = P * S;

    % Define constant orbital parameters for all satellites
    eccentricity = zeros(1, total_sats);
    argumentOfPerigee = zeros(1, total_sats);
    semiMajorAxis = repmat(semiMajorAxis, 1, total_sats);
    inclination = repmat(inclination, 1, total_sats);
    RAAN = repelem((0:P-1) * 180 / P, S);
    trueAnomaly = repmat((0:S-1) * 360 / S, 1, P);

    % Create the satellite array
    sat = satellite(sc, semiMajorAxis, eccentricity, inclination, RAAN, argumentOfPerigee, trueAnomaly);

    % Define the Ground Station (GS)
    gs = groundStation(sc, gs_lat, gs_long);
    gs_ecef_pos = lla2ecef([gs_lat, gs_long, gs_alt]);
    gs_ecef_pos_time = repmat(gs_ecef_pos', 1, num_times);

    % Define dimensions
    num_sats = numel(sat);

    % Define array to store positions of satellites through time (3D: 3 coords × times × satellites)
    all_sat_pos = states(sat, "CoordinateFrame", "ecef");

    % Define array for Access status to signal visibility between GS and sats
    ac = access(sat, gs);
    intvls = accessIntervals(ac);

    % Initialize status array to store logical access status at each time
    status = false(num_times, num_sats);

    % Extract status for each time step
    for j = 1:num_times
        status(j, :) = accessStatus(ac, times(j));
    end

    % Define the number of closest satellites to find
    n = 3;

    % Preallocate matrices to store data for n closest satellites
    closest_sat_coords = zeros(num_times, n, 3);
    closest_sat_indices = zeros(num_times, n);
    closest_sat_dists = zeros(num_times, n);
    closest_sat_elevations = zeros(num_times, n);
    closest_sat_elevations_discrete = zeros(num_times, n);

    % Find the closest satellites and calculate distances
    fprintf('Calculating closest satellites and distances...\n');
    for j = 1:num_times
        active_sats = find(status(j, :));
        gs_pos = gs_ecef_pos_time(:, j);

        if ~isempty(active_sats)
            distances = zeros(1, length(active_sats));
            sat_positions_active = zeros(3, length(active_sats));

            for k = 1:length(active_sats)
                sat_idx = active_sats(k);
                sat_pos = all_sat_pos(:, j, sat_idx);
                distances(k) = norm(gs_pos - sat_pos);
                sat_positions_active(:, k) = sat_pos;
            end

            [sorted_dists, dist_idx] = sort(distances);
            num_closest = min(n, length(active_sats));
            closest_sats_indices_in_active = dist_idx(1:num_closest);
            closest_sats = active_sats(closest_sats_indices_in_active);
            closest_dists = sorted_dists(1:num_closest);

            for k = 1:num_closest
                sat_idx = closest_sats(k);
                sat_pos = sat_positions_active(:, closest_sats_indices_in_active(k));
                closest_sat_coords(j, k, :) = sat_pos;
                closest_sat_indices(j, k) = sat_idx;
                closest_sat_dists(j, k) = closest_dists(k);
            end
        end

        if mod(j, round(num_times/10)) == 0 || j == num_times
            fprintf('Processed %d/%d time steps (%.0f%%)\n', j, num_times, (j/num_times)*100);
        end
    end

    % Compute and store elevation angles for the closest satellites
    fprintf('Calculating elevation angles for closest satellites...\n');
    R_E = 6371000; % Earth's mean radius in meters
    h_0 = semiMajorAxis(1) - R_E;
    gs_lat_rad = deg2rad(gs_lat);
    gs_long_rad = deg2rad(gs_long);

    for j = 1:num_times
        closest_sats = closest_sat_indices(j, :);
        num_closest = sum(closest_sats > 0);

        for k = 1:num_closest
            sat_idx = closest_sats(k);
            if sat_idx > 0
                sat_pos = squeeze(closest_sat_coords(j, k, :));
                ssp_lla = ecef2lla(sat_pos');
                ssp_lla(3) = 0;
                ssp_lat_rad = deg2rad(ssp_lla(1));
                ssp_long_rad = deg2rad(ssp_lla(2));
                delta_phi = gs_long_rad - ssp_long_rad;
                gamma_ell = acos(sin(ssp_lat_rad) * sin(gs_lat_rad) + cos(ssp_lat_rad) * cos(gs_lat_rad) * cos(delta_phi));

                if sin(gamma_ell) > eps
                    term = (cos(gamma_ell) - (R_E / (R_E + h_0))) / sin(gamma_ell);
                    elevation = rad2deg(atan(term));
                else
                    elevation = 90;
                end

                closest_sat_elevations(j, k) = elevation;
            end
        end

        if mod(j, round(num_times/10)) == 0 || j == num_times
            fprintf('Elevation processed %d/%d time steps (%.0f%%)\n', j, num_times, (j/num_times)*100);
        end
    end

    % Discretize elevation angles
    fprintf('Discretizing elevation angles...\n');
    min_elevation = 15;
    discretized_values = [20 30 45 60 70];
    
    for j = 1:num_times
        for k = 1:n
            elevation = closest_sat_elevations(j, k);

            if elevation > min_elevation
                [~, closest_idx] = min(abs(discretized_values - elevation));
                closest_sat_elevations_discrete(j, k) = discretized_values(closest_idx);
            else
                closest_sat_elevations_discrete(j, k) = 0;
            end
        end
    end

    fprintf('Discretization complete.\n');
    gs.MarkerColor = [1 0 0]; % Red marker for ground station
    play(sc);
end
