%% simple pulse

clc;clear;close all; 
% Define the parameters of the rectangular pulse
pulse_duration = 1; % Pulse duration in seconds
Ts = 0.05;          % Time resolution in seconds


% Create the time and Doppler shift vectors
time = -2*pulse_duration:Ts:2*pulse_duration;

% Generate the rectangular pulse
rect_pulse = rectpuls(time, pulse_duration);

% Call my amiguity function
[t_pulse,f_pulse,ambg_pulse] = ambgfunn(rect_pulse,time);


%% LFM Signal
clc;close all;
% Define parameters
Fs = 1000; % Sampling frequency
T = 1; % Duration of the LFM signal in seconds
t = 0:1/Fs:T-1/Fs; % Time vector

% Define different slopes for LFM signals
slopes = [-100, -50, 50, 100]; % in Hz/second



% Generate and plot LFM signals with different slopes
figure;
for i = 1:numel(slopes)
    slope = slopes(i);
    signal = exp(1j * pi * slope * t.^2);
    
    subplot(numel(slopes), 1, i);
    [t_lfm,f_lfm,ambg_lfm] = ambgfunn(signal, t);
    title(sprintf('Ambiguity Function of LFM Signal (Slope: %d Hz/second)', slope));
end




%% Barker code
% Generate Barker code signal
barker_code5 = [1, 1, 1, -1, 1];
barker_code7 = [1 1 1 -1 -1 1 -1];
barker_code11 =  [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
barker_code13 = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];

% Chage here to chage output
barker_code = barker_code13;

% Define time vector
time_vector = linspace(0, 1, length(barker_code));

% Call the ambgfunn function
[t_barker,f_barker,ambg_barker] = ambgfunn(barker_code, time_vector);




%% Compare sidelobe levels 
ambg_pulse_zeroFrequency = ambg_pulse(floor(length(f_pulse)),:);
ambg_barker_zeroFrequency = ambg_barker(floor(length(f_barker)),:);
ambg_lfm_zeroFrequency = ambg_lfm(floor(length(f_lfm)),:);


figure
hold on
plot(t_lfm,mag2db(ambg_lfm_zeroFrequency))
grid on
xlabel('Delay \tau (second)')
ylabel('Ambiguity Function (dB)')
title('AF Barker Code, 0 Hz Doppler Cut')



%% 5
% Define parameters
Fs = 1000; % Sampling frequency
T = 1; % Duration of the signals in seconds
t = 0:1/Fs:T-1/Fs; % Time vector

% Generate LFM signals with positive and negative slopes
slope_positive = 50; % Positive slope in Hz/second
slope_negative = -50; % Negative slope in Hz/second

lfm_signal_positive = exp(1j * pi * slope_positive * t.^2);
lfm_signal_negative = exp(1j * pi * slope_negative * t.^2);

% Plot the ambiguity functions of the LFM signals
figure;
subplot(2, 1, 1);
ambgfunn(lfm_signal_positive, t);
title('Ambiguity Function of LFM Signal (Positive Slope)');

subplot(2, 1, 2);
ambgfunn(lfm_signal_negative, t);
title('Ambiguity Function of LFM Signal (Negative Slope)');



%% 
function [t, f, ambiguity] = ambgfunn(signal,time_vector)
    %This is my ambiguity function
    %   signal is pulse and time_vector is time vector including times

    L = length(signal);
    Ts = time_vector(2)-time_vector(1); % Get time sampling
    f = 1/Ts*(   floor(-L/2):floor(L/2)-1)/L; % Define frequency vector
    
    % Initialize output
    ambiguity = zeros(L, 2*length(time_vector));
    % Delay the signal 2*its length (+-)
    for i = -L:L
        s1 = signal;
        % s2 ~~ s*(t-taw)
        s2 = circshift(signal,i);
        if i>0
            s2(1:i) = 0;
        elseif i<0
            s2(end+i+1:end) = 0;
        end
        s2 = conj(s2);
        % Multiply in time domain
        s = s1.*s2;
        % Get FFT
        sf = fft(s);
        sf = fftshift(sf);
        df = flip(sf);
        % Put the vector for specific taw in matrix
        ambiguity(:,i+L+1) = abs(df);
    end

    % Plot the ambiguity function matrix
    contour((-L:L)*Ts,f, ambiguity);
    xlabel('Time (second)');
    ylabel('Doppler Shift (hertz)');
    title('Ambiguity Function of Rectangular Pulse');
    colorbar;  

    t = (-L:L)*Ts;

end

