function [Sat_Ang_time, full_fadedWave, full_inputWave, full_channel, full_state, full_time] = WP3(Sample_Rate, num_times, closest_sat_elevations_discrete)

    seed = 72;

    all_faded_waves = cell(num_times, 1);
    all_channel_gains = cell(num_times, 1);
    all_states = cell(num_times, 1);
    Sat_Ang_time = zeros(num_times, 1);
    
    new_State = "good";

    for t = 1:num_times
        elev = closest_sat_elevations_discrete(t);
        Sat_Ang_time(t) = elev;

        if elev > 0
            chan = p681LMSChannel;
            chan.Environment = "Urban";
            chan.CarrierFrequency = 3.8e9;
            chan.ElevationAngle = elev;
            chan.MobileSpeed = 0.5;
            chan.AzimuthOrientation = 0;
            chan.SampleRate = Sample_Rate;
            chan.InitialState = new_State;
            chan.FadingTechnique = "Filtered Gaussian noise";
            chan.Seed = seed + t;

            rng(seed);

            chanDur = 60;
            in = complex(randn(numSamples, 1), randn(numSamples, 1));
            [fadedWave, channelCoefficients, sampleTimes, stateSeries] = step(chan, in);

            all_faded_waves{t} = fadedWave;
            all_channel_gains{t} = channelCoefficients;
            all_states{t} = stateSeries;
            all_times{t} = sampleTimes;

            previousState = stateSeries(end);
            new_State = "good" if previousState ~= 0 else "bad";
        else
            all_faded_waves{t} = NaN;
            all_channel_gains{t} = NaN;
            all_states{t} = NaN;
            all_times{t} = NaN;
        end
    end

    full_time = [];
    full_fadedWave = [];
    full_inputWave = [];
    full_channel = [];
    full_state = [];
    t_offset = 0;

    for idx = max(1, num_times - 20 + 1):num_times
        fadedWave = all_faded_waves{idx};
        channelCoefficients = all_channel_gains{idx};
        stateSeries = all_states{idx};

        full_time = [full_time; sampleTimes + t_offset];
        full_fadedWave = [full_fadedWave; fadedWave];
        full_channel = [full_channel; channelCoefficients];
        full_state = [full_state; stateSeries];

        t_offset = t_offset + max(sampleTimes);
    end

    Last_20_Ang = Sat_Ang_time(max(1, end-19:end));
    fprintf('Last 20 Elevation Angles: ');
    fprintf('%f ', Last_20_Ang);
    fprintf('\n');

    figure(1)
    plot(full_time, 20*log10(abs(full_inputWave)), full_time, 20*log10(abs(full_fadedWave)))
    title('Power Profile - Last 20 Simulations')
    legend('Input Waveform', 'Faded Waveform')
    xlabel('Time (s)')
    ylabel('Power (dB)')

    figure(2)
    plot(full_time, 20*log10(abs(full_channel)))
    title('Channel Gain - Last 20 Simulations')
    xlabel('Time (s)')
    ylabel('Path Gain (dB)')

    figure(3)
    plot(full_time, full_state)
    title('State Series - Last 20 Simulations')
    axis([0 full_time(end) -0.5 1.5])
    xlabel('Time (s)')
    ylabel('State')

    figure(4)
    plot(Sat_Ang_time, 'o-', 'MarkerSize', 3, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
    title('Satellite Elevation Angle over Time')
    xlabel('Time Step')
    ylabel('Elevation Angle (degrees)')

    figure(5)
    histogram(Sat_Ang_time(~isnan(Sat_Ang_time)), 'BinWidth', 5, 'FaceColor', 'g')
    title('Elevation angles occurrency')
    xlabel('Elevation (degrees)')
    ylabel('Occurencies')
    grid on
end
