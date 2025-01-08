clear all
close all

Fs = 8192; %sampling frequency

notecreate = @(notenumber,dur) sin(2*pi*[1:floor(dur*Fs)]/Fs * 440*2.^((notenumber-1)/12));

%piano notes in sequence starting from A4 (440 Hz)
notename = {'A' 'A#' 'B' 'C' 'C#' 'D' 'D#' 'E' 'F' 'F#' 'G' 'G#'}; 

% simple song
song = {'A' 'A' 'E' 'E' 'F#' 'F#' 'E' 'E' 'D' 'D' 'C#' 'C#' 'B' 'B' 'A' 'A'};

% find piano note number 
songidx = zeros(1,length(song));
for k1 = 1:length(song)
    idx = strcmp(song(k1), notename);
    songidx(k1) = find(idx);
end    

dur = 0.3;% Duration of each note in seconds

% write notes into the song
songnotes = [];
for k1 = 1:length(songidx)
    y_note = notecreate(songidx(k1),dur);
    y_note = [y_note, zeros(1,75)]';
    songnotes = [songnotes; y_note];
    
end

% plot the song
t = 0:1/Fs:(length(songnotes)-1)/Fs; %time
figure, plot(t,songnotes)

% play song
soundsc(songnotes, Fs)