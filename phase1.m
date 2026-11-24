clear; close all; clc;

cfg = config();

if ~exist(cfg.PHASE1_DIR, 'dir')
    mkdir(cfg.PHASE1_DIR);
end

exts  = {'*.wav', '*.mp3', '*.m4a', '*.flac', '*.ogg'};
files = [];
for k = 1:numel(exts)
    files = [files; dir(fullfile(cfg.INPUT_DIR, exts{k}))]; 
end
if isempty(files)
    error('No sound files found in folder "%s".', cfg.INPUT_DIR);
end

plotFile = cfg.PLOT_FILE;
if isempty(plotFile)
    plotFile = files(1).name;
end
[~, plotBase, ~] = fileparts(plotFile);

sounds = struct('name', {}, 'x', {}, 'fs', {});
for k = 1:numel(files)
    fname      = fullfile(cfg.INPUT_DIR, files(k).name);
    [~, baseName, ~] = fileparts(files(k).name);
    isPlotFile = strcmpi(baseName, plotBase);
    doWave     = cfg.PLOT_ALL || isPlotFile;   % Task 3.5 waveform plot
    % doCos    = isPlotFile;                   % Task 3.7 cosine plot (commented out for Phase 2)

    [x, fs] = preprocessSound(fname, cfg.PHASE1_DIR, cfg.TARGET_FS, ...
                              cfg.PLAY_AUDIO, doWave);

    sounds(k).name = files(k).name;
    sounds(k).x    = x;
    sounds(k).fs   = fs;
end
fprintf('\nPhase 1 complete: %d file(s) saved to "%s".\n', numel(files), cfg.PHASE1_DIR);

function [x, fs] = preprocessSound(fname, outDir, targetFs, playAudio, doWave)

[~, name, ~] = fileparts(fname);
fprintf('\n=== %s ===\n', name);

% ---- Task 3.1: read the file and find its sampling rate ----
[x, fsOrig] = audioread(fname);
fprintf('Sampling rate: %d Hz\n', fsOrig);

% ---- Task 3.2: stereo or mono? ----
[nSamples, nChannels] = size(x);
fprintf('Channels: %d | Samples: %d | Duration: %.2f s\n', ...
        nChannels, nSamples, nSamples / fsOrig);
if nSamples / fsOrig > 30
    warning('%s is longer than 30 s (the handout asks for 30 s or less).', name);
end
if nChannels == 2
    x = x(:, 1) + x(:, 2);           % add the two columns -> one column
    fprintf('Stereo: the two channels were added to make a mono signal.\n');
else
    fprintf('Mono: no conversion needed.\n');
end

% ---- Task 3.6: resample to 16 kHz ----
if fsOrig < targetFs
    warning(['%s was recorded at %d Hz, below %d Hz. Upsampling cannot add ' ...
             'the missing high frequencies, so re-record it at a higher rate.'], ...
            name, fsOrig, targetFs);
end
if fsOrig ~= targetFs
    x = resampleFFT(x, targetFs, fsOrig);
    fprintf('Resampled from %d Hz to %d Hz (%d samples).\n', fsOrig, targetFs, length(x));
end
fs = targetFs;

if max(abs(x)) > 1
    x = x / max(abs(x));
end

% ---- Task 3.3: play the sound ----
if playAudio
    fprintf('Playing %s ...\n', name);
    playblocking(audioplayer(x, fs));
end

% ---- Task 3.4: write the sound to a new file ----
outFile = fullfile(outDir, [name '_mono16k.wav']);
audiowrite(outFile, x, fs);
fprintf('Written: %s\n', outFile);

% ---- Task 3.5: waveform vs sample number ----
if doWave
    fig = figure('Name', ['Task 3.5 - ' name]);
    plot(1:length(x), x);
    xlabel('Sample number');
    ylabel('Amplitude');
    title(sprintf('Waveform of %s (mono, %d Hz)', name, fs), 'Interpreter', 'none');
    xlim([1 length(x)]);
    grid on;
    print(fig, fullfile(outDir, ['task3_5_waveform_' name '.png']), '-dpng', '-r200');
end

% ---- Task 3.7: Commented out for Phase 2 ----
% cosFreq = 1000;
% N       = length(x);
% t       = (0:N-1)' / fs;
% cosSig  = cos(2*pi*cosFreq*t);
% fprintf('Cosine: %d Hz, %d samples, %.2f s\n', cosFreq, N, N / fs);
%
% if doCos
%     if playAudio
%         fprintf('Playing the %d Hz cosine ...\n', cosFreq);
%         playblocking(audioplayer(0.3 * cosSig, fs));
%     end
%
%     T    = 1 / cosFreq;
%     nTwo = round(2 * T * fs) + 1;
%     fig  = figure('Name', 'Task 3.7 - cosine');
%     plot(t(1:nTwo) * 1000, cosSig(1:nTwo), '-o', 'MarkerSize', 4);
%     xlabel('Time (ms)');
%     ylabel('Amplitude');
%     title(sprintf('Two cycles of the %d Hz cosine (f_s = %d Hz)', cosFreq, fs));
%     xlim([0 2 * T * 1000]);
%     ylim([-1.1 1.1]);
%     grid on;
%     print(fig, fullfile(outDir, 'task3_7_cosine.png'), '-dpng', '-r200');
% end
end

function y = resampleFFT(x, fsOut, fsIn)
% RESAMPLEFFT  Resample a column vector from fsIn to fsOut using the FFT.
%   Keeping only frequencies below the new Nyquist limit acts as an ideal
%   anti-aliasing low-pass filter.
x = x(:);
N = length(x);
M = round(N * fsOut / fsIn);
X = fft(x);
Y = zeros(M, 1);
h = floor((min(N, M) - 1) / 2);
Y(1:h+1)       = X(1:h+1);            % DC and positive frequencies
Y(end-h+1:end) = X(end-h+1:end);      % negative frequencies
y = real(ifft(Y)) * (M / N);
end
