%% clear workspace
clc
clear all
close all
addpath 'matlab-midi/src' % use midi library
%% Data Import
% loading time series (can be anything: seismic, potential fields, climate record)
X = load('data');
% X = load('W01-20-ascii.dat');

dt = 4e-3; % sample  interval (s)
%dt = 1/20;
Fs = 1/dt; % sampling frequency (Hz)
f_noise = 50/(Fs/2);

% X(:,5)=[];
% X(1:400000,:)=[];
% %X=X(1:end/3,:);
% maxx = max(abs(X));
% X(:,1)=(X(:,1)-median(X(:,1)))/maxx(1);
% X(:,2)=(X(:,2)-median(X(:,2)))/maxx(2);
% X(:,3)=(X(:,3)-median(X(:,3)))/maxx(3);
% X(:,4)=(X(:,4)-median(X(:,4)))/maxx(4);

% Design a Notch Filter to remove the noise at f_noise
% wo = f_noise/(Fs/2);   % Normalized Frequency
% bw = wo/35;            % Bandwidth for the notch
% [b, a] = iirnotch(wo, bw); % Notch filter design
% f1=1/100/(Fs/2);
% f2=1/3/(Fs/2);
f1=3/(Fs/2)
f2=60/(Fs/2)
% Design a Bandpass Filter to isolate frequencies between f1 and f2
[b, a] = butter(4, [f1, f2]); % 4th order Butterworth bandpass filter

% Apply Filter
X(:,1) = filter(b, a, X(:,1));
% X(:,4) = filter(b, a, X(:,4));

X = X(:,1);
% X(:,1)=2*X(:,1);
nt = length(X);
t = (0:(nt-1))*dt;
%%
figure, plot(t,X), xlabel('Time(s)'),ylabel('Amplitude'), axis tight
%%
% writematrix([t',X(:,1)],'datatest.csv')
%% convert to WAV (Waveform Audio File Format)
% audio file format for storing waveform data. Essentially, it's a digital representation of sound.
% raw, uncompressed audio data. This data represent the actual sound waves.
% WAV files are large, detailed information about the sound wave.
% used for high-quality audio recordings like music tracks, sound effects, and voice recordings.
% widely compatible with different audio playback and editing software.
% File name for the output WAV file
%y = resample(50*X,t,500*Fs);
%figure, plot(t,X), xlabel('Time(s)'),ylabel('Amplitude'), axis tight
%hold on, plot(0:1/100/Fs:1/100/Fs*(length(y)-1),y,'r')
filename = 'data.wav';
% Write the signal to a WAV file
% audiowrite(filename, X, 44100);
%%
y = audioread(filename);
% Fs1 = 41000;
% dt1 = 1/8000;
% t1 = 0:dt1:max(t);
% y = interp1(t, y, t1); 
sound(y, 44100);

%% Compute spectrogram using short-time Fourier transform
% Parameters for the STFT
windowSize = 256;      % Size of the window
overlap = round(0.75 * windowSize); % Overlap between successive windows (75% overlap)
nfft = 256;            % Number of FFT points
% Compute the spectrogram
[S, F, T] = stft(X, Fs, 'Window', hamming(windowSize), ...
    'OverlapLength', overlap, 'FFTLength', nfft);
% spectrum is symmetric, we take a half
 F = F(nfft/2:nfft);
 S =  S(nfft/2:nfft,:);
% Convert to power 
P = abs(S).^2;
% Convert power to decibels
S_dB = 10*log10(abs(P/1e-12));

% Display the spectrogram
figure;
subplot(211),plot(t,X),axis tight
subplot(212),surf(T, F, S_dB, 'EdgeColor', 'none');
axis xy; axis tight; colormap(jet); view(0, 90);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of the Signal');
colorbar('southoutside'); ylim([0 120])
%%
k = 3;% number of maximum values
S_dB = medfilt2(S_dB,[3 3]); % get amplitude spectrum and smooth spectrogram using median filter
% Find peak frequencies in each time slice
[maxVal,Ind]=sort(S_dB,'descend');
maxVal = maxVal(1:k,:);
Ind = Ind(1:k,:);
peakFrequencies = zeros(size(T,1),k);
for i = 1:k
 peakFrequencies(:,i) = F(Ind(i,:));% find peak frequency
end 
FreqCoeff = 200; % scaling frequency coefficient to get audible signal  
%
figure, imagesc(T,F,abs(S_dB)), hold on, 
for i = 1:k
plot(T,peakFrequencies(:,i),'+')
end
ylim([0 120])
%% Convert frequencies to MIDI note numbers
midiNotes = ceil(58 + 12 * log2(FreqCoeff*peakFrequencies / 440));
% initialize matrix:
N = size(midiNotes,1);  % number of notes
%loudness ('velocity') of signal normalized in range 0 to 127 (2^7)
loudness = 127*max(abs(S)/max(abs(S(:))));
%% Create MIDI structure and write into disk
k=1;
M = zeros(N*k,1);
for i = 1:k
    ii = (i-1)*N+(1:N);
    M(ii,1) = 1;         % track 1
    M(ii,2) = i;         % channel 1
    M(ii,3) = midiNotes(:,i)+(i-1)*3; % note numbers: one octave starting at middle C (60)
    M(ii,4) = loudness;
    M(ii,5) = 0.3*(1:N);  % note on:  notes start every .3 seconds
    M(ii,6) = M(ii,5) + .3;   % note off: each note has duration .3 seconds
end

midi_new = matrix2midi(M); % convert matrix to MIDI structure
writemidi(midi_new, 'AudioSeismogram-1.mid'); % write MIDI file onn disk
%% Plot and play sonified seismic trace
figure('WindowState','maximized')% create and maximize figure
y = midi2audio(midi_new,Fs*FreqCoeff,'fm');  % convert MIDI structure to audio to play in MATLAB
player = audioplayer(y, Fs*FreqCoeff); % create audio object
%play(player); %start playing 
for it = 1:nt
    subplot(311) % 3-rows and 1-column panel figure
    plot(t,X)% plot signal in the upper plane 
    hold on, 
    plot([T(it), T(it)],[min(X) max(X)],'k','LineWidth',2) % plot pointer for time 
    hold off
    ylabel('Normalize trace amplitude (n.d.)')
    subplot(312)
    imagesc(T,F,S_dB), colormap(flipud(bone))% instaneous spectrum plot
    hold on,
    plot([T(it), T(it)],[min(F) max(F)],'k','LineWidth',2)
    hold off
    ylim([0 50])
    ylabel('Frequency (Hz)')
    subplot(313)
    plot(T,loudness(:)), hold on,
    plot([T(it), T(it)],[0 150],'k','LineWidth',2)
    hold off
    xlabel('Trace time (s)'), ylabel('Sound amplitude (pnt)')
    drawnow
    pause(0.01)
end
%
stop(player)