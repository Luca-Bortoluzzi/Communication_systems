%% Exercise 3 - Part 3
clc;
clear;
close all;
%Values from previous part 2 to make the code run
B = 36000000;
SNR_dB = 20.7152;

%% DVB-S2X MODCOD Lookup Table (ETSI EN 302 307-2, Slide 40)
modcod_names = { ...
    'QPSK 2/9'; 'QPSK 13/45'; 'QPSK 9/20'; 'QPSK 11/20'; ...
    '8APSK 5/9-L'; '8APSK 26/45-L'; '8PSK 23/36'; '8PSK 25/36'; '8PSK 13/18'; ...
    '16APSK 1/2-L'; '16APSK 8/15-L'; '16APSK 5/9-L'; '16APSK 26/45'; ...
    '16APSK 3/5'; '16APSK 3/5-L'; '16APSK 28/45'; '16APSK 23/36'; ...
    '16APSK 2/3-L'; '16APSK 25/36'; '16APSK 13/18'; '16APSK 7/9'; '16APSK 77/90'; ...
    '32APSK 2/3-L'; '32APSK 32/45'; '32APSK 11/15'; '32APSK 7/9'; ...
    '64APSK 32/45-L'; '64APSK 11/15'; '64APSK 7/9'; '64APSK 4/5'; '64APSK 5/6'; ...
    '128APSK 3/4'; '128APSK 7/9'; ...
    '256APSK 29/45-L'; '256APSK 2/3-L'; '256APSK 31/45-L'; ...
    '256APSK 32/45'; '256APSK 11/15-L'; '256APSK 3/4'};

eta_spec_vec = [ ...
    0.434841; 0.567805; 0.889135; 1.088581; ...
    1.647211; 1.713601; 1.896173; 2.062148; 2.145136; ...
    1.972253; 2.104850; 2.193247; 2.281645; ...
    2.370043; 2.370043; 2.458441; 2.524739; ...
    2.635236; 2.745734; 2.856231; 3.077225; 3.386618; ...
    3.291954; 3.510192; 3.620536; 3.841226; ...
    4.206428; 4.338659; 4.603122; 4.735354; 4.936639; ...
    5.163248; 5.355556; ...
    5.065690; 5.241514; 5.417338; ...
    5.593162; 5.768987; 5.900855];

req_EsN0_dB = [ ...
    -2.85; -2.03; 0.22; 1.45; ...
    4.73;  5.13; 6.12; 7.02; 7.49; ...
    5.97;  6.55; 6.84; 7.51; ...
    7.80;  7.41; 8.10; 8.38; ...
    8.43;  9.27; 9.71; 10.65; 11.99; ...
    11.10; 11.75; 12.17; 13.05; ...
    13.98; 14.81; 15.47; 15.87; 16.55; ...
    17.73; 18.53; ...
    16.98; 17.24; 18.10; ...
    18.59; 18.84; 19.57];

%% 3.3.1
% Under alpha = 0: Rs = B and Es/N0 = SNR
EsN0_avail_dB = SNR_dB;

% Filter all MODCODs whose required Es/N0 does not exceed the available Es/N0
valid_idx = find(req_EsN0_dB <= EsN0_avail_dB);

if isempty(valid_idx)
    error('Available Es/N0 is below the minimum supported MODCOD threshold.');
end

% Select the supported MODCOD with the maximum spectral efficiency
[~, best_sub_idx] = max(eta_spec_vec(valid_idx));
best_idx = valid_idx(best_sub_idx);

selected_MODCOD = modcod_names{best_idx};
eta_best        = eta_spec_vec(best_idx);
req_EsN0_best   = req_EsN0_dB(best_idx);
link_margin_dB  = EsN0_avail_dB - req_EsN0_best;

%% 3.3.2-3-4
Rs = B;                                 % Symbol rate [sym/s] (alpha = 0)
Rb = Rs * eta_best;                     % Achievable bit rate [bit/s]

R_HD = 5e6;                             % Average bit rate per HD channel [bit/s]
N_HD_channels = floor(Rb / R_HD);       % Integer number of simultaneous HD channels

SNR_lin = 10^(SNR_dB / 10);             % Linear SNR
C = B * log2(1 + SNR_lin);              % Shannon channel capacity [bit/s]

efficiency_vs_shannon = (Rb / C) * 100; % Percentage of Shannon limit achieved

%% 3.3.5
figure;

% Bar chart comparing Achievable Bit Rate (Rb) and Shannon Capacity (C)
b = bar([Rb / 1e6, C / 1e6], 0.5);
b.FaceColor = 'flat';
b.CData(1, :) = [0.00, 0.45, 0.74]; % Blue for Rb
b.CData(2, :) = [0.85, 0.33, 0.10]; % Orange-red for C

grid on;
set(gca, 'XTickLabel', {sprintf('Achievable Bit Rate R_b\n(%s)', selected_MODCOD), ...
    'Shannon Capacity C'});
ylabel('Data Rate [Mbit/s]');
title(sprintf('DVB-S2X Bit Rate vs. Shannon Limit (B = %.0f MHz, SNR = %.2f dB)', ...
    B / 1e6, SNR_dB));
ylim([0, (C / 1e6) * 1.18]);

% Add numerical labels above each bar
text(1, Rb / 1e6 + 3, sprintf('%.2f Mbit/s\n(%.1f%% of C | %d HD ch.)', ...
    Rb / 1e6, efficiency_vs_shannon, N_HD_channels), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(2, C / 1e6 + 3, sprintf('%.2f Mbit/s\n(Theoretical Bound)', C / 1e6), ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% Comments about the obtained results
%    Proximity to the Shannon Limit (~84.3% Efficiency):
%    With an available SNR (= Es/N0) of 13.52 dB over B = 36 MHz, the highest
%    supported DVB-S2X MODCOD is 32APSK 7/9 (required Es/N0 = 13.05 dB,
%    spectral efficiency = 3.841226 bit/symbol), yielding Rb = 138.28 Mbit/s.
%    Compared to the theoretical Shannon capacity C = 163.95 Mbit/s, the
%    system achieves 84.34% of the channel limit, confirming the strong error-
%    correcting performance of the concatenated LDPC + BCH codes in DVB-S2X.
%
%    MODCOD Discretization and Residual Link Margin:
%    Unlike Shannon capacity C, which scales continuously with SNR, a practical
%    receiver must lock onto a discrete MODCOD profile. Because our SNR
%    (13.52 dB) falls 0.46 dB short of the next entry in the table (64APSK
%    32/45-L, requiring 13.98 dB), the system operates with a 0.47 dB link
%    margin above the 32APSK 7/9 threshold. While this margin leaves some
%    spectral capacity unused under clear-sky conditions, it provides valuable
%    resilience against minor atmospheric or pointing fluctuations.
%
%    Idealized Modeling Assumptions and Broadcast Capacity:
%    At 5 Mbit/s per channel, Rb = 138.28 Mbit/s supports 27 simultaneous HD
%    television channels (with 3.28 Mbit/s of residual capacity). In an actual
%    operational link, however, a non-zero roll-off factor (typically alpha =
%    0.05 to 0.20), physical-layer framing overhead (PLFRAME headers and phase-
%    tracking pilots needed for 32APSK), and transport stream encapsulation
%    would slightly reduce the net user throughput.