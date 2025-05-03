function sat = WP2_create_sat_constellation(sc, P, S, semiMajorAxis, inclination)
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