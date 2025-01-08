%% 
clear all
close all
% Parameters
Fs = 44100; % Sampling frequency in Hz
T = 5;      % Duration of the signal in seconds
f0 = 500;   % Start frequency of the chirp (500 Hz)
f1 = 3000;  % End frequency of the chirp (3000 Hz)

% Time vector
t = 0:1/Fs:T;

% Generate the chirp signal
x = (1-t/T).^2.*chirp(t, f0, T, f1);
% Play the chirp signal
sound(x, Fs);
% Plot the chirp signal
plot(t, x);
xlabel('Time (s)');
ylabel('Amplitude');
title('Chirp Signal');
axis([0 T -1.1 1.1]); % Set axis limits for better visualization



%%
% x = load('data2');
% dt = 4e-3; % sample  interval (s)
%Fs = 1/dt; % sampling frequency (Hz)
dt = 1/Fs;
nt = length(x);
%t = (0:(nt-1))*dt;
% [x,t] = resample(x,t,Fs*10);
% Parameters for the STFT
windowSize = 2^11;      % Size of the window
overlap = round(0.75 * windowSize); % Overlap between successive windows (75% overlap)
nfft = 2^11;            % Number of FFT points
% Compute the spectrogram
[S, F, T] = stft(x, Fs, 'Window', hamming(windowSize), ...
    'OverlapLength', overlap, 'FFTLength', nfft);
% Convert to power 
P = abs(S).^2;
% Convert power to decibels
S_dB = 10*log10(abs(P/1e-12));
% Display the spectrogram
figure;
subplot(211),plot(t,x),axis tight
subplot(212),surf(T, F(nfft/2:nfft), S_dB(nfft/2:nfft,:), 'EdgeColor', 'none');
axis xy; axis tight; colormap(jet); view(0, 90);ylim([0 5000])
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of the Signal');
yy=colorbar('southoutside');
ylabel(yy,'dB')
%%
% Parameters
Fs = 44100; % Sampling frequency in Hz
T = 5;      % Duration of the signal in seconds
f0 = 500;   % Start frequency of the chirp (500 Hz)
f1 = 3000;  % End frequency of the chirp (3000 Hz)

% Time vector
t = 0:1/Fs:T;

% Generate the chirp signal
x = (1-t/T).^2.*chirp(t, f0, T, f1);
% Perform Continuous Wavelet Transform
[wt, f] = cwt(x, Fs,'FrequencyLimits',[0 5000],'VoicesPerOctave',10);

% Plot the Spectrogram using Wavelet Transform
figure;
surface(t, f, abs(wt));
axis tight;
shading flat;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram using Continuous Wavelet Transform');
colorbar;
