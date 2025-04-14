%%

% Parameters
seed = 72;
Sample_Rate = 400;
numSamples = Sample_Rate * num_times * 60; % Sample rate times duration of simulation in seconds

all_faded_waves = cell(num_times, 1);
all_channel_gains = cell(num_times, 1);
all_states = cell(num_times, 1);
all_distances = cell(num_times, 1);

Sat_Ang_time = zeros(num_times,1);

% Loop over time steps
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
        chan.InitialState = "Good";
        chan.FadingTechnique = "Filtered Gaussian noise";
        chan.RandomStream = "mt19937ar with seed";
        chan.Seed = seed + t; % ensure unique seed for each time

        % Set random number generator with seed
        rng(seed);
        % Channel duration 60 sec, because we update the angle every 60 sec
        chanDur = 60; 
        % Random input waveform
        numSamples = floor(chan.SampleRate*chanDur)+1;
        in = complex(randn(numSamples,1),randn(numSamples,1));
        % Pass the input signal through channel
        [fadedWave,channelCoefficients,sampleTimes,stateSeries] = step(chan,in);

        % Store
        all_faded_waves{t} = fadedWave;
        all_channel_gains{t} = channelCoefficients;
        all_states{t} = stateSeries;
        all_times{t} = sampleTimes;
    else
        % No satellite in view — store NaN
        all_faded_waves{t} = NaN;
        all_channel_gains{t} = NaN;
        all_states{t} = NaN;
        all_times{t} = NaN;
    end
end



% Find the last valid simulation index
valid_idx = find(~cellfun(@(x) any(isnan(x(:))), all_faded_waves), 1, 'last');


fadedWave = all_faded_waves{valid_idx};
channelCoefficients = all_channel_gains{valid_idx};
stateSeries = all_states{valid_idx};
distanceVector = all_times{valid_idx};

timeVector = all_times{valid_idx};  

figure(1)
plot(timeVector, 20*log10(abs(in(1:length(timeVector)))), timeVector, 20*log10(abs(fadedWave)))
title(['Power Profile of Waveform over Time (Speed = ' num2str(chan.MobileSpeed) ' m/s)'])
legend('Input Waveform', 'Faded Waveform')
xlabel('Time (in seconds)')
ylabel('Power (in dB)')

figure(2)
plot(timeVector, 20*log10(abs(channelCoefficients)))
title('Channel Gain over Time')
xlabel('Time (in seconds)')
ylabel('Path Gain (in dB)')

figure(3)
plot(timeVector, stateSeries)
title('State Series of Channel over Time')
axis([0 timeVector(end) -0.5 1.5])
xlabel('Time (in seconds)')
ylabel('State')

figure(4)
plot(Sat_Ang_time)
title('Satellite Elevation Angle over Time')
xlabel('Time Step')
ylabel('Elevation Angle (degrees)')
