function = WP3(Sample_Rate, numSamples)

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

end
