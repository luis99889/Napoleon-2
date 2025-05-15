
function [Bit_rate, channel_gain_avg, L] = WP4_function_cooperative(all_channel_gains, closest_sat_dists, closest_sat_indices, CarrierFrequency, cooperative_mode)
% WP4_function_cooperative: Calculate bit rate with satellite cooperation
% Inputs:
%   all_channel_gains: Cell array of channel gains
%   closest_sat_dists: Distances to closest satellites
%   closest_sat_indices: Indices of closest satellites (needed for cooperation)
%   CarrierFrequency: Carrier frequency
%   cooperative_mode: 0=no cooperation, 1=pair cooperation, 2=full cooperation

% Constants
f = CarrierFrequency;
c = 3e8;
L = zeros(1, length(closest_sat_dists));
G_gs = 1;
G_sat = 10^5;
N0 = 3.98e-21;
BW = 5e6;
N = length(closest_sat_dists);
P_Tx = 10e-3; % 10 mW

% Free-space path loss
for i = 1:length(closest_sat_dists)
    L(i) = ((c)/(4*pi*closest_sat_dists(i,1)*f))^2; 
end

% Channel gain averaging
channel_gain_avg = zeros(1, length(all_channel_gains));
for j = 1:length(all_channel_gains)
    channel_gain_avg(j) = mean(abs(all_channel_gains{j}))^2;
end

% Calculate received power for each satellite
PRX = zeros(1, length(closest_sat_dists));
for z = 1:length(closest_sat_dists) 
    PRX(z) = P_Tx * G_sat * G_gs * channel_gain_avg(z) * L(z);
end

% Calculate bit rate based on cooperative mode
Bit_rate = 0;

if cooperative_mode == 0
    % Non-cooperative mode (original implementation)
    for k = 1:length(all_channel_gains)
        Bit_rate = Bit_rate + (BW * (log2(1 + PRX(k)/(N0 * BW))));
    end
    Bit_rate = Bit_rate / N;
    
elseif cooperative_mode == 1
    % Pair-wise cooperation (use best 2 satellites at each time)
    % Assuming time dimension is handled elsewhere and we're processing a single time step
    
    % Sort satellites by received power
    [sorted_PRX, sort_idx] = sort(PRX, 'descend');
    
    % If we have at least 2 satellites
    if length(PRX) >= 2
        % Combine power from the best two satellites
        combined_PRX = sorted_PRX(1) + sorted_PRX(2);
        
        % Calculate bit rate with combined SNR
        Bit_rate = BW * log2(1 + combined_PRX/(N0 * BW));
        
        fprintf('Cooperating satellites: %d and %d\n', ...
                closest_sat_indices(sort_idx(1)), closest_sat_indices(sort_idx(2)));
    else
        % If only one satellite available
        Bit_rate = BW * log2(1 + PRX(1)/(N0 * BW));
    end
    
elseif cooperative_mode == 2
    % Full cooperation (all satellites)
    % Sum received power from all available satellites
    combined_PRX = sum(PRX);
    
    % Calculate bit rate with combined SNR
    Bit_rate = BW * log2(1 + combined_PRX/(N0 * BW));
    
    fprintf('Using all %d satellites in cooperation\n', length(PRX));
end

% Display the calculated bit rate
disp(['Bit rate: ', num2str(Bit_rate), ' bits/s']);

end

