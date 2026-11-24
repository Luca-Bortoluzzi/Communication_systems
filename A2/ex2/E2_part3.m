%% Exercise 2 - Part 3: CDMA Downlink (2 users, no cyclic shift)
clc;
clear;
close all;

%% 1. System Parameters
N = 12;                         % Spreading factor (sequence length)
Eb = N;                         % Energy per bit (N unit-energy chips)
target_errors = 100;            % Minimum number of bit errors per simulated SNR point

EbN0_sim_db = 0:8;              % Simulated Eb/N0 range [dB] (Point 11)
EbN0_sim_lin = 10.^(EbN0_sim_db / 10);

EbN0_th_db = 0:0.5:20;          % Theoretical Eb/N0 range [dB]
EbN0_th_lin = 10.^(EbN0_th_db / 10);

% Analytical 2-PAM benchmark curve (Point 10)
P_2pam = 0.5 * erfc(sqrt(EbN0_th_lin));

rng(42);                        % Fixed random seed for reproducibility

% Generate random binary sequence c1 for U1 (same as in Part 2)
c1 = randi([0, 1], 1, N);
c1_bipolar = 2*c1 - 1;

%% 2. Loop over the 4 required values of parameter p
for p = [0, 2, 4, 12]

    % Build c2 by flipping (N - p)/2 bits of c1 so that |dot(c1, c2)| == p
    c2 = c1;
    num_flips = (N - p) / 2;
    c2(1:num_flips) = 1 - c2(1:num_flips);
    c2_bipolar = 2*c2 - 1;

    % Print sequences to Command Window to include them in the report (Point 10)
    fprintf('p = %2d | c1 = [%s] | c2 = [%s]\n', p, num2str(c1), num2str(c2));

    % --- Theoretical curves for the current value of p ---
    % Exact analytical curve (Point 12)
    P_analytic = 0.25 * erfc(sqrt(EbN0_th_lin) * (1 + p/N)) + ...
        0.25 * erfc(sqrt(EbN0_th_lin) * (1 - p/N));

    % High-SNR approximate curve accounting for multiplicity (Point 13)
    if p == 0
        P_approx = 0.5 * erfc(sqrt(EbN0_th_lin));
    else
        P_approx = 0.25 * erfc(sqrt(EbN0_th_lin) * (1 - p/N));
    end

    % --- Monte Carlo simulation from 0 to 8 dB (Point 11) ---
    BER_sim = zeros(size(EbN0_sim_db));

    for i = 1:length(EbN0_sim_db)
        errors = 0;
        bits = 0;

        N0 = Eb / EbN0_sim_lin(i);
        sigma = sqrt(N0 / 2);

        while errors < target_errors
            % Generate random bits and convert to bipolar (-1, +1)
            b1 = randi([0, 1]);
            b1_bip = 2*b1 - 1;

            b2 = randi([0, 1]);
            b2_bip = 2*b2 - 1;

            % Spread (encode) both users' bits and sum with AWGN noise
            s1 = b1_bip * c1_bipolar;
            s2 = b2_bip * c2_bipolar;
            noise = sigma * randn(1, N);

            r = s1 + s2 + noise;

            % U1 receiver: de-spread using c1_bipolar and make decision
            y = dot(r, c1_bipolar);
            a1 = (y >= 0);          % 1 if y >= 0, 0 if y < 0

            if a1 ~= b1
                errors = errors + 1;
            end
            bits = bits + 1;
        end

        BER_sim(i) = errors / bits;
    end

    % --- Plot all 4 curves for the current p (Points 10-13) ---
    figure;
    theme(gcf, "light");
    semilogy(EbN0_th_db, P_2pam, 'r--', 'LineWidth', 2); hold on;
    semilogy(EbN0_th_db, P_analytic, 'b-', 'LineWidth', 1.5);
    semilogy(EbN0_sim_db, BER_sim, 'ko-', 'LineWidth', 1.2, 'MarkerSize', 7);
    semilogy(EbN0_th_db, P_approx, 'm-.', 'LineWidth', 1.5);

    ylim([1e-8 1]);
    xlim([0 20]);
    grid on;

    set(gca, 'FontSize', 12);
    xlabel('E_b/N_0 [dB]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('BER', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('CDMA Downlink, 2 users, N = %d, p = %d', N, p), 'FontSize', 14, 'FontWeight', 'bold');
    legend('2-PAM', 'Analytic', 'Simulation', 'Approximated', 'Location', 'northeast', 'FontSize', 11);
end

%% Point 14 - Comments on the Simulation & Analytical Curves
% - Across all tested values of p (0, 2, 4, 12), the Monte Carlo simulation
%   points (0 to 8 dB) match the exact analytical BER curve, validating the model.
% - For p = 0, the exact and approximate curves coincide over the entire Eb/N0
%   range and overlap with the single-user 2-PAM benchmark.
% - For p = 2 and p = 4, the high-SNR approximate curve slightly underestimates
%   the exact BER at low Eb/N0 (0 to 4 dB) because it neglects constructive
%   interference, but converges to the exact curve at higher Eb/N0 (>= 6 dB).
% - For p = 12, both the exact and approximate curves flatten out into a
%   horizontal line at BER = 0.25.

%% Point 15 - Impact of Parameter p on CDMA Downlink Performance
% - The parameter p = |c1' * c2'| quantifies the Multiple Access Interference (MAI)
%   from user U2 leaking into U1's receiver. Increasing p progressively degrades BER.
% - When U1 and U2 transmit opposite bits (50% probability), destructive
%   interference reduces the useful signal amplitude from N to (N - p), causing
%   an asymptotic SNR penalty of -20*log10(1 - p/N) dB:
%     * p = 0  (orthogonal): 0 dB loss -> MAI is cancelled, achieving 2-PAM BER.
%     * p = 2  (low MAI):    1.58 dB SNR penalty at high Eb/N0.
%     * p = 4  (moderate):   3.52 dB SNR penalty at high Eb/N0.
%     * p = 12 (identical):  total signal cancellation (1 - p/N = 0) when bits
%                            differ, creating an irreducible error floor at
%                            BER = 0.25 regardless of how high Eb/N0 is.
% - System takeaway: since downlink signals are synchronous (zero relative delay),
%   assigning mutually orthogonal sequences (p = 0, e.g., Walsh-Hadamard codes)
%   completely eliminates intra-cell interference without extra transmit power.