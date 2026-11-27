%% Exercise 2 - Part 4: CDMA Uplink (2 users, with cyclic shifts)
clc;
clear;
close all;

%% Point 16 - Asymptotic Governing Parameter p (CDMA Uplink)
% - In the asynchronous CDMA uplink, user U2's signal arrives with a random
%   cyclic shift tau uniformly distributed in {0, 1, ..., N-1} at each bit.
% - At high Eb/N0 (asymptotic regime), the average BER is dominated by the
%   worst-case cyclic shift producing the strongest destructive interference.
% - Therefore, the parameter p governing asymptotic performance is the maximum
%   cyclic cross-correlation magnitude (r_C) across all N cyclic shifts:
%       p = max_{0 <= tau <= N-1} |c1' · T^tau(c2')|

%% Point 17 - Exact Analytical BER Curve (with Cyclic Shifts)
% - Received vector: r = b1' * c1' + b2' * T^tau(c2') + n, where b1', b2' in {-1, +1}
%   are equiprobable and tau is uniformly distributed in {0, ..., N-1} (prob. 1/N).
% - Despreading at U1 receiver yields:
%       y = r · c1' = b1' * N + b2' * r_{c1',c2'}(tau) + w,   w ~ N(0, N * N0/2)
% - Defining p_tau = |r_{c1',c2'}(tau)| and averaging over both equiprobable
%   bits b2' and all N equiprobable shifts tau (2N combinations, prob. 1/(2N)):
%       P_b(e) = (1 / (4*N)) * sum_{tau=0}^{N-1} [ ...
%           erfc(sqrt(Eb/N0) * (1 + p_tau / N)) + ...
%           erfc(sqrt(Eb/N0) * (1 - p_tau / N)) ]

%% Point 18 - High-SNR Approximate BER Curve
% - At high Eb/N0, the sum is dominated by the worst-case destructive interference
%   term (1 - p/N), where p = max_tau(p_tau).
% - Let Kp be the number of cyclic shifts tau in {0, ..., N-1} achieving p_tau = p.
%   Out of 2N equiprobable (b2', tau) pairs, the worst-case event occurs Kp times
%   (multiplicity Kp / (2N)):
%       * For p > 0: P_approx(e) = (Kp / (4*N)) * erfc(sqrt(Eb/N0) * (1 - p/N))
%       * For p = 0: P_approx(e) = 0.5 * erfc(sqrt(Eb/N0))  (orthogonal for all tau)

%% 1. System Parameters
N = 12;                         % Spreading factor
Eb = N;                         % Energy per bit
target_errors = 100;            % Minimum bit errors per simulated SNR point

EbN0_sim_db = 0:8;              % Simulated Eb/N0 range [dB]
EbN0_sim_lin = 10.^(EbN0_sim_db / 10);

EbN0_th_db = 0:0.5:20;          % Theoretical Eb/N0 range [dB]
EbN0_th_lin = 10.^(EbN0_th_db / 10);

% Analytical 2-PAM benchmark curve (Point 19)
P_2pam = 0.5 * erfc(sqrt(EbN0_th_lin));

rng(42);                        % Fixed seed for reproducibility

% Generate balanced c1 for U1 (guarantees existence of a partner even for p = 0)
c1 = [0 0 0 0 0 0 1 1 1 1 1 1];
c1_bipolar = 2*c1 - 1;

%% --- plot 2PAM ---
figure;
theme(gcf, "light");
semilogy(EbN0_th_db, P_2pam, 'r--', 'LineWidth', 2);
ylim([1e-8 1]); xlim([0 20]); grid on;
set(gca, 'FontSize', 12);
xlabel('E_b/N_0 [dB]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('BER', 'FontSize', 12, 'FontWeight', 'bold');
title('2-PAM Analytical Benchmark', 'FontSize', 14, 'FontWeight', 'bold');
legend('2-PAM', 'Location', 'northeast', 'FontSize', 11);
exportgraphics(gcf, "A2_2_4_2PAM.png", Resolution=300);

for target_p = [0, 2, 4, 12]

    % Draw random binary c2 (without balance constraint, required for p = 2)
    while true
        c2 = randi([0, 1], 1, N);
        c2_bipolar = 2*c2 - 1;

        % Periodic cross-correlation via FFT
        R_c1c2 = round(real(ifft(fft(c1_bipolar) .* conj(fft(c2_bipolar)))));
        R_c1c2 = fftshift(R_c1c2);
        p_tau = abs(R_c1c2);

        if max(p_tau) == target_p
            break;
        end
    end

    p = max(p_tau);             % Asymptotic governing parameter (Point 16)
    Kp = sum(p_tau == p);       % Multiplicity of worst-case shift (Point 18)

    fprintf('p = %2d (Kp = %2d) | c1 = [%s] | c2 = [%s]\n', ...
        p, Kp, num2str(c1), num2str(c2));

    % Precompute all N cyclic shifts of c2_bipolar to speed up Monte Carlo
    C2_shifts = zeros(N, N);
    for tau = 0:N-1
        C2_shifts(tau + 1, :) = circshift(c2_bipolar, tau);
    end

    % --- Exact analytical curve averaged over all N cyclic shifts (Point 21) ---
    P_analytic = zeros(size(EbN0_th_lin));
    for k = 1:N
        P_analytic = P_analytic + (1/N) * ( ...
            0.25 * erfc(sqrt(EbN0_th_lin) * (1 + p_tau(k)/N)) + ...
            0.25 * erfc(sqrt(EbN0_th_lin) * (1 - p_tau(k)/N)) );
    end

    % --- High-SNR approximate curve accounting for multiplicity Kp/(2N) (Point 22) ---
    if p == 0
        P_approx = 0.5 * erfc(sqrt(EbN0_th_lin));
    else
        P_approx = (Kp / (4*N)) * erfc(sqrt(EbN0_th_lin) * (1 - p/N));
    end

    % --- Fast Monte Carlo simulation from 0 to 8 dB (Point 20) ---
    BER_sim = zeros(size(EbN0_sim_db));
    c1_col = c1_bipolar';       % Column vector for fast matrix multiplication

    for i = 1:length(EbN0_sim_db)
        errors = 0;
        bits = 0;

        N0 = Eb / EbN0_sim_lin(i);
        sigma = sqrt(N0 / 2);

        while errors < target_errors
            b1 = randi([0, 1]);
            b1_bip = 2*b1 - 1;

            b2_bip = 2*randi([0, 1]) - 1;
            idx_shift = randi(N);   % Random cyclic shift index in 1..N

            r = b1_bip * c1_bipolar + b2_bip * C2_shifts(idx_shift, :) + sigma * randn(1, N);

            y = r * c1_col;
            a1 = (y >= 0);

            if a1 ~= b1
                errors = errors + 1;
            end
            bits = bits + 1;
        end

        BER_sim(i) = errors / bits;
    end
    
    %% 3. Progressive Plots (Points 19 to 22)
    
    % % --- Plot 2 (Point 20): 2-PAM + Simulation ---
    % figure;
    % theme(gcf, "light");
    % semilogy(EbN0_th_db, P_2pam, 'r--', 'LineWidth', 2); hold on;
    % semilogy(EbN0_sim_db, BER_sim, 'ko', 'LineWidth', 1.2, 'MarkerSize', 7);
    % ylim([1e-8 1]); xlim([0 20]); grid on;
    % set(gca, 'FontSize', 12);
    % xlabel('E_b/N_0 [dB]', 'FontSize', 12, 'FontWeight', 'bold');
    % ylabel('BER', 'FontSize', 12, 'FontWeight', 'bold');
    % title(sprintf('CDMA uplink, 2 users, N = %d, p = %d (Step 2: + Simulation)', N, p), 'FontSize', 13, 'FontWeight', 'bold');
    % legend('2-PAM', 'simulation', 'Location', 'northeast', 'FontSize', 11);
    % 
    % % --- Plot 3 (Point 21): 2-PAM + Simulation + Analytic ---
    % figure;
    % theme(gcf, "light");
    % semilogy(EbN0_th_db, P_2pam, 'r--', 'LineWidth', 2); hold on;
    % semilogy(EbN0_th_db, P_analytic, 'b-', 'LineWidth', 1.5);
    % semilogy(EbN0_sim_db, BER_sim, 'ko', 'LineWidth', 1.2, 'MarkerSize', 7);
    % ylim([1e-8 1]); xlim([0 20]); grid on;
    % set(gca, 'FontSize', 12);
    % xlabel('E_b/N_0 [dB]', 'FontSize', 12, 'FontWeight', 'bold');
    % ylabel('BER', 'FontSize', 12, 'FontWeight', 'bold');
    % title(sprintf('CDMA uplink, 2 users, N = %d, p = %d (Step 3: + Analytic)', N, p), 'FontSize', 13, 'FontWeight', 'bold');
    % legend('2-PAM', 'analytic', 'simulation', 'Location', 'northeast', 'FontSize', 11);
    
    % --- Plot 4 (Point 22): All 4 curves ---
    figure;
    theme(gcf, "light");
    semilogy(EbN0_th_db, P_2pam, 'r--', 'LineWidth', 2); hold on;
    semilogy(EbN0_th_db, P_analytic, 'b-', 'LineWidth', 1.5);
    semilogy(EbN0_sim_db, BER_sim, 'ko', 'LineWidth', 1.2, 'MarkerSize', 7);
    semilogy(EbN0_th_db, P_approx, 'r--', 'LineWidth', 1);
    ylim([1e-8 1]); xlim([0 20]); grid on;
    set(gca, 'FontSize', 12);
    xlabel('E_b/N_0 [dB]', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('BER', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('CDMA uplink, 2 users, N = %d, p = %d', N, p), 'FontSize', 13, 'FontWeight', 'bold');
    legend('2-PAM', 'analytic', 'simulation', 'approximated', 'Location', 'northeast', 'FontSize', 11);
    fileName = sprintf("A2_2_4_%d.png", target_p);
    exportgraphics(gcf, fileName, "resolution", 300);
end