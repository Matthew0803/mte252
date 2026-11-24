clear; close all; clc;

cfg = config();

if ~exist(cfg.PHASE1_DIR, 'dir')
    error('Phase 1 folder "%s" not found. Run phase1.m first.', cfg.PHASE1_DIR);
end
if ~exist(cfg.PHASE2_DIR, 'dir')
    mkdir(cfg.PHASE2_DIR);
end

files = dir(fullfile(cfg.PHASE1_DIR, '*_mono16k.wav'));
if isempty(files)
    error(['No *_mono16k.wav files found in "%s". ' ...
           'Run phase1.m first to generate Phase 1 outputs.'], cfg.PHASE1_DIR);
end

% Resolve which file to plot if cfg.PLOT_ALL == false
plotBase = stripMonoSuffix(cfg.PLOT_FILE);
if isempty(plotBase)
    plotBase = stripMonoSuffix(files(1).name);
end

useSPT = hasSignalToolbox();
if useSPT
    fprintf('Filter backend: Signal Processing Toolbox (butter/zp2sos/sosfilt)\n');
else
    fprintf('Filter backend: Base MATLAB bilinear-transform SOS (no toolbox needed)\n');
end

bpBank = designBandpassBank(cfg.N_CHANNELS, cfg.F_MIN, cfg.F_MAX, ...
                            cfg.TARGET_FS, cfg.BP_ORDER, cfg.SPACING, useSPT);
lpf    = designEnvelopeLPF(cfg.LPF_CUTOFF, cfg.TARGET_FS, cfg.LPF_ORDER, useSPT);

% Save frequency-response plots to document filter designs (Tasks 4 & 8)
plotFilterResponses(bpBank, lpf, cfg.TARGET_FS, cfg.PHASE2_DIR, useSPT);

sounds = struct('name', {}, 'x', {}, 'fs', {}, 'bpOut', {}, 'env', {});
for k = 1:numel(files)
    fname    = fullfile(cfg.PHASE1_DIR, files(k).name);
    baseName = stripMonoSuffix(files(k).name);
    fprintf('\n=== %s ===\n', baseName);

    % Load preprocessed Phase 1 audio directly
    [x, fs] = audioread(fname);
    if size(x, 2) ~= 1
        error('%s is not mono (%d channels). Re-run phase1.m.', files(k).name, size(x, 2));
    end
    if fs ~= cfg.TARGET_FS
        error('%s has fs = %d Hz (expected %d Hz). Re-run phase1.m.', ...
              files(k).name, fs, cfg.TARGET_FS);
    end
    fprintf('Loaded Phase 1 file: %s (%d samples, %.2f s, %d Hz)\n', ...
            files(k).name, length(x), length(x) / fs, fs);

    doPlots = cfg.PLOT_ALL || strcmpi(baseName, plotBase);

    % Tasks 5 - 9: bandpass filtering, rectification, LPF envelope, plots
    [bpOut, env] = processPhase2(x, fs, baseName, cfg.PHASE2_DIR, bpBank, lpf, ...
                                 doPlots, useSPT);

    sounds(k).name  = baseName;
    sounds(k).x     = x;
    sounds(k).fs    = fs;
    sounds(k).bpOut = bpOut;
    sounds(k).env   = env;
end

if cfg.SAVE_MAT
    matFile = fullfile(cfg.PHASE2_DIR, 'phase2_results.mat');
    save(matFile, 'sounds', 'bpBank', 'lpf', 'cfg', '-v7.3');
    fprintf('\nSaved Phase 2 workspace data to "%s".\n', matFile);
end

fprintf('\nPhase 2 complete: %d file(s) processed. Outputs in "%s".\n', ...
        numel(files), cfg.PHASE2_DIR);


% Phase 2 functions
function base = stripMonoSuffix(filename)
% STRIPMONOSUFFIX  Removes extension and trailing '_mono16k' for clean naming.
if isempty(filename)
    base = '';
    return;
end
[~, base, ~] = fileparts(filename);
suffix = '_mono16k';
if length(base) > length(suffix) && strcmpi(base(end-length(suffix)+1:end), suffix)
    base = base(1:end-length(suffix));
end
end

function bpBank = designBandpassBank(nChannels, fMin, fMax, fs, order, spacing, useSPT)
% Task 4: Design an N-channel Butterworth bandpass bank.
% Calls designBandpassFilter() inside a for loop for each channel.

if strcmpi(spacing, 'log')
    edges = logspace(log10(fMin), log10(fMax), nChannels + 1);
elseif strcmpi(spacing, 'linear')
    edges = linspace(fMin, fMax, nChannels + 1);
else
    error('cfg.SPACING must be ''log'' or ''linear''.');
end

bpBank.nChannels = nChannels;
bpBank.spacing   = spacing;
bpBank.order     = order;
bpBank.edges     = edges;
bpBank.fCenter   = zeros(1, nChannels);
bpBank.sos       = cell(1, nChannels);
bpBank.g         = zeros(1, nChannels);

fprintf('\n--- Task 4: Bandpass Filter Bank (%d channels, %s spacing, order %d) ---\n', ...
        nChannels, spacing, order);
for ch = 1:nChannels
    f1 = edges(ch);
    f2 = edges(ch + 1);

    if strcmpi(spacing, 'log')
        bpBank.fCenter(ch) = sqrt(f1 * f2);   % geometric center (for Phase 3)
    else
        bpBank.fCenter(ch) = (f1 + f2) / 2;   % arithmetic center (for Phase 3)
    end

    [bpBank.sos{ch}, bpBank.g(ch)] = designBandpassFilter(f1, f2, fs, order, useSPT);
    fprintf('  Ch %2d: %7.1f Hz - %7.1f Hz  (fc = %7.1f Hz)\n', ...
            ch, f1, f2, bpBank.fCenter(ch));
end
end

function [sos, g] = designBandpassFilter(fLow, fHigh, fs, order, useSPT)
% Single IIR Butterworth bandpass filter in SOS form.
% Normalized frequency Wn = f / (fs/2) must lie strictly inside (0, 1).

nyquist = fs / 2;
if fLow <= 0 || fHigh >= nyquist || fLow >= fHigh
    error('Band edges [%.1f, %.1f] Hz must satisfy 0 < fLow < fHigh < %.1f Hz.', ...
          fLow, fHigh, nyquist);
end

Wn = [fLow, fHigh] / nyquist;         % normalized frequency (pi rad/sample)

if useSPT
    [z, p, k] = butter(order, Wn, 'bandpass');
    [sos, g]  = zp2sos(z, p, k);
else
    [sos, g]  = butterBandSOSManual(order, Wn(1), Wn(2));
end
end

function lpf = designEnvelopeLPF(fCutoff, fs, order, useSPT)
% Task 8: Butterworth lowpass filter for envelope detection.

nyquist = fs / 2;
if fCutoff <= 0 || fCutoff >= nyquist
    error('LPF cutoff %.1f Hz must lie strictly inside (0, %.1f) Hz.', fCutoff, nyquist);
end

Wn = fCutoff / nyquist;               % normalized cutoff frequency

if useSPT
    [z, p, k] = butter(order, Wn, 'low');
    [lpf.sos, lpf.g] = zp2sos(z, p, k);
else
    [lpf.sos, lpf.g] = butterLowSOSManual(order, Wn);
end

lpf.fCutoff = fCutoff;
lpf.order   = order;
fprintf('--- Task 8: Envelope LPF (Butterworth, f_c = %d Hz, order %d) ---\n', ...
        fCutoff, order);
end

function [bpOut, env] = processPhase2(x, fs, name, outDir, bpBank, lpf, doPlots, useSPT)
% Runs Tasks 5 - 9 for one preprocessed 16 kHz mono signal.

nSamples  = length(x);
nChannels = bpBank.nChannels;
t         = (0:nSamples-1)' / fs;

% ---- Task 5: Filter the sound with the bandpass filter bank ----
bpOut = zeros(nSamples, nChannels);
for ch = 1:nChannels
    bpOut(:, ch) = applySOS(bpBank.sos{ch}, bpBank.g(ch), x, useSPT);
end

% ---- Task 6: Plot outputs of the lowest and highest frequency channels ----
if doPlots
    fig = figure('Name', ['Task 6 - BPF Outputs (' name ')'], ...
                 'Position', [100 100 900 560]);

    subplot(2, 1, 1);
    plot(t, bpOut(:, 1));
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Lowest channel (Ch 1: %.1f - %.1f Hz)', ...
                  bpBank.edges(1), bpBank.edges(2)));
    xlim([0 t(end)]);
    grid on;

    subplot(2, 1, 2);
    plot(t, bpOut(:, nChannels));
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Highest channel (Ch %d: %.1f - %.1f Hz)', ...
                  nChannels, bpBank.edges(nChannels), bpBank.edges(nChannels + 1)));
    xlim([0 t(end)]);
    grid on;

    sgtitle(sprintf('Task 6: Bandpass Filter Outputs - %s', name), ...
            'Interpreter', 'none');
    outImg = fullfile(outDir, ['task6_bpf_outputs_' name '.png']);
    print(fig, outImg, '-dpng', '-r200');
    fprintf('Saved plot: %s\n', outImg);
end

% ---- Task 7: Envelope extraction step 1 (full-wave rectification) ----
rectOut = abs(bpOut);

% ---- Task 8: Envelope extraction step 2 (400 Hz lowpass filter) ----
% Implemented manually (no built-in envelope() function).
env = zeros(nSamples, nChannels);
for ch = 1:nChannels
    y = applySOS(lpf.sos, lpf.g, rectOut(:, ch), useSPT);
    env(:, ch) = max(y, 0);   % clamp negligible negative ringing near zero
end

% ---- Task 9: Plot extracted envelopes of the lowest and highest channels ----
if doPlots
    fig = figure('Name', ['Task 9 - Envelopes (' name ')'], ...
                 'Position', [140 140 900 560]);

    subplot(2, 1, 1);
    plot(t, env(:, 1), 'LineWidth', 1.1);
    xlabel('Time (s)');
    ylabel('Envelope amplitude');
    title(sprintf('Lowest channel envelope (Ch 1: %.1f - %.1f Hz)', ...
                  bpBank.edges(1), bpBank.edges(2)));
    xlim([0 t(end)]);
    grid on;

    subplot(2, 1, 2);
    plot(t, env(:, nChannels), 'LineWidth', 1.1);
    xlabel('Time (s)');
    ylabel('Envelope amplitude');
    title(sprintf('Highest channel envelope (Ch %d: %.1f - %.1f Hz)', ...
                  nChannels, bpBank.edges(nChannels), bpBank.edges(nChannels + 1)));
    xlim([0 t(end)]);
    grid on;

    sgtitle(sprintf('Task 9: Extracted Envelopes (LPF %d Hz) - %s', ...
                    lpf.fCutoff, name), 'Interpreter', 'none');
    outImg = fullfile(outDir, ['task9_envelopes_' name '.png']);
    print(fig, outImg, '-dpng', '-r200');
    fprintf('Saved plot: %s\n', outImg);
end
end

function plotFilterResponses(bpBank, lpf, fs, outDir, useSPT)
% PLOTFILTERRESPONSES  Saves magnitude responses for the Phase 2 report.

nPts = 4096;

% Task 4: Bandpass bank magnitude response
fig1 = figure('Name', 'Task 4 - Bandpass Filter Bank Response', ...
              'Position', [80 80 900 420]);
hold on;
for ch = 1:bpBank.nChannels
    [H, fAxis] = evalSOSFreqResp(bpBank.sos{ch}, bpBank.g(ch), nPts, fs, useSPT);
    plot(fAxis, 20 * log10(abs(H) + eps), 'LineWidth', 1.2);
end
hold off;
grid on;
xlim([0 fs/2]);
ylim([-60 5]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title(sprintf('Task 4: %d-Channel Butterworth Bandpass Bank (%s spacing, order %d)', ...
              bpBank.nChannels, bpBank.spacing, bpBank.order));
print(fig1, fullfile(outDir, 'task4_filterbank_response.png'), '-dpng', '-r200');

% Task 8: Envelope LPF magnitude response
fig2 = figure('Name', 'Task 8 - Envelope LPF Response', ...
              'Position', [120 120 900 420]);
[Hlpf, fAxis] = evalSOSFreqResp(lpf.sos, lpf.g, nPts, fs, useSPT);
plot(fAxis, 20 * log10(abs(Hlpf) + eps), 'LineWidth', 1.3);
grid on;
xlim([0 2000]);
ylim([-60 5]);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title(sprintf('Task 8: Butterworth Envelope Lowpass Filter (f_c = %d Hz, order %d)', ...
              lpf.fCutoff, lpf.order));
print(fig2, fullfile(outDir, 'task8_lpf_response.png'), '-dpng', '-r200');
end


function tf = hasSignalToolbox()
% Returns true only if Signal Processing Toolbox is licensed and available.
tf = license('test', 'Signal_Toolbox') && ...
     exist('butter', 'file') > 0 && ...
     exist('zp2sos', 'file') > 0 && ...
     exist('sosfilt', 'file') > 0 && ...
     exist('freqz', 'file') > 0;
end

function y = applySOS(sos, g, x, useSPT)
% Applies a cascaded second-order section filter with overall gain g.
if useSPT
    y = g * sosfilt(sos, x);
else
    y = x;
    for s = 1:size(sos, 1)
        y = filter(sos(s, 1:3), sos(s, 4:6), y);
    end
    y = g * y;
end
end

function [H, fAxis] = evalSOSFreqResp(sos, g, nPts, fs, useSPT)
% Evaluates the complex frequency response H(f) on [0, fs/2].
if useSPT
    [H, fAxis] = freqz(sos, nPts, fs);
    H = g * H;
else
    fAxis = linspace(0, fs/2, nPts)';
    zInv  = exp(-1j * 2 * pi * fAxis / fs);
    H     = ones(nPts, 1) * g;
    for s = 1:size(sos, 1)
        num = sos(s, 1) + sos(s, 2) * zInv + sos(s, 3) * (zInv .^ 2);
        den = sos(s, 4) + sos(s, 5) * zInv + sos(s, 6) * (zInv .^ 2);
        H   = H .* (num ./ den);
    end
end
end

function [sos, g] = butterLowSOSManual(order, Wn)
% Exact bilinear-transform Butterworth lowpass filter in SOS form.
K   = tan(pi * Wn / 2);
sos = [];

for k = 1:floor(order / 2)
    theta = pi * (2*k + order - 1) / (2*order);
    s     = K * exp(1j * theta);
    zp    = (1 + s) / (1 - s);
    sos   = [sos; 1, 2, 1, 1, -2*real(zp), abs(zp)^2]; 
end
if mod(order, 2) == 1
    s   = -K;
    zp  = (1 + s) / (1 - s);
    sos = [sos; 1, 1, 0, 1, -real(zp), 0]; 
end

% Normalize DC gain (z = 1) to unity
H0 = 1;
for s = 1:size(sos, 1)
    H0 = H0 * sum(sos(s, 1:3)) / sum(sos(s, 4:6));
end
g = real(1 / H0);
end

function [sos, g] = butterBandSOSManual(order, W1, W2)
% Exact bilinear-transform Butterworth bandpass filter in SOS form.
Om1    = tan(pi * W1 / 2);
Om2    = tan(pi * W2 / 2);
B      = Om2 - Om1;
Om0_sq = Om1 * Om2;
sos    = [];

for k = 1:floor(order / 2)
    theta  = pi * (2*k + order - 1) / (2*order);
    pk     = exp(1j * theta);
    disc   = sqrt((B * pk)^2 - 4 * Om0_sq);
    sPair  = [(B * pk + disc) / 2, (B * pk - disc) / 2];
    for idx = 1:2
        zp  = (1 + sPair(idx)) / (1 - sPair(idx));
        sos = [sos; 1, 0, -1, 1, -2*real(zp), abs(zp)^2]; 
    end
end
if mod(order, 2) == 1
    pk   = -1;
    disc = sqrt((B * pk)^2 - 4 * Om0_sq);
    s1   = (B * pk + disc) / 2;
    s2   = (B * pk - disc) / 2;
    zp1  = (1 + s1) / (1 - s1);
    zp2  = (1 + s2) / (1 - s2);
    sos  = [sos; 1, 0, -1, 1, -real(zp1 + zp2), real(zp1 * zp2)]; 
end

% Normalize gain at warped center frequency to unity
w0   = 2 * atan(sqrt(Om0_sq));
zInv = exp(-1j * w0);
H0   = 1;
for s = 1:size(sos, 1)
    num = sos(s, 1) + sos(s, 2) * zInv + sos(s, 3) * (zInv ^ 2);
    den = sos(s, 4) + sos(s, 5) * zInv + sos(s, 6) * (zInv ^ 2);
    H0  = H0 * (num / den);
end
g = 1 / abs(H0);
end
