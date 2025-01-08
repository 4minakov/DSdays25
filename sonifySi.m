% Sonification of Silicon Electron Configuration Including 3s and 3p Orbitals
clear all
close all
% Binding energies in eV (estimated for 3s and 3p)
energy_1s = 1839;   % K (1s)
energy_2s = 149.7;  % L1 (2s)
energy_2p = (99.82 + 99.42) / 2; % Average for L2 (2p1/2) and L3 (2p3/2)
energy_3s = 50;     % Estimated for 3s
energy_3p = 45;     % Estimated for 3p

% Scaling factor (arbitrary for illustration)
scaling_factor = 2; 

% Convert energies to frequencies within audible range
freq_1s = energy_1s * scaling_factor;
freq_2s = energy_2s * scaling_factor;
freq_2p = energy_2p * scaling_factor;
freq_3s = energy_3s * scaling_factor;
freq_3p = energy_3p * scaling_factor;

% Electron configuration of Silicon: 
e_conf = {'1s²', '2s²', '2p⁶', '3s²', '3p²'};
electrons_in_orbitals = [2, 2, 6, 2, 2]; % 1s, 2s, 2p, 3s, 3p
frequencies = [freq_1s, freq_2s, freq_2p, freq_3s, freq_3p];

% Duration of each note
duration = 0.5; % 0.5 seconds

% Sampling frequency
fs = 44.1e3; % Standard sampling rate
fs = 8192; % Standard sampling rate

% Generate and play the sounds
for i = 1:length(electrons_in_orbitals)
    % Generate a sine wave for the current frequency
    t = 0:1/fs:duration;
    y = sin(2 * pi * frequencies(i) * t) + 0.5*sin(2 * pi * 2 * frequencies(i) * t)...
        + 0.25*sin(2 * pi * 3 * frequencies(i) * t) + 0.125*sin(2 * pi * 4 * frequencies(i) * t);
    figure(1),clf
    subplot(211)
    plot(t,y),title(e_conf{i}),xlabel('time'),ylabel('amplitude')
    subplot(212)
    [wt, f] = cwt(y,fs,'FrequencyLimits',[0 1000]);
    S_dB = 10*log10(abs(abs(wt).^2/1e-12));
    surface(t, f, abs(wt).^2);
    axis tight;
    shading flat;
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Spectrogram using Continuous Wavelet Transform');
%     colorbar
    drawnow
    % Play the sound for the number of electrons in the orbital
    for j = 1:electrons_in_orbitals(i)
        sound(y, fs);
        pause(duration + 0.3);
    end 
end


% End of script


