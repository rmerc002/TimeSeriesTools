function [d sr duration] = readAudio(usedSampleRate)
%% Open selected audion from all the formats .mp3 .bwf .wav .ogg .aiff
% Convert it to wav format
% Downsamle it to SampleRate
% Convert stereo to mono

[filename, pathname] = uigetfile( {'*.mp3;*.bwf;*.wav;*.ogg;*.AIFF;*All Files (*.*)'},'Pick a file');
fullPath = strcat(pathname,filename)
[ y, Fs ] = audioread( fullPath );
if(strfind(filename, '.mp3')) % Convert to Wav Format
    wav_file = 'convertWavFile.wav';
    audiowrite(wav_file,y,Fs);
else
    wav_file = fullPath;
end

info = audioinfo(wav_file);
duration = getfield(info,'Duration');
samplerate = getfield(info,'SampleRate');

[ d, sr ] = audioread( wav_file );
[m , n] = size(d); % convert stereo to monot
if(n > 1)
    d = d(:,1)+d(:,2);
end
if(samplerate ~= usedSampleRate) % change sample rate
    [P,Q] = rat(samplerate/usedSampleRate);
    d = resample(d,P,Q); 
end

end