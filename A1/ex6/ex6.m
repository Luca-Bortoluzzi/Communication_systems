clc;
close all;
clear;
%% --- Properties of m sequences --- 
% Design a randomizer with m = 10
m = 10; % number of flip-flops
% chosen polynomial
% P_x = x^{10} + x^7 + 1;
% choose a random seed
state = [1 1 1 1 1 0 0 0 0 0];
% generate an entire m-sequence
N = 2^m - 1;
prbs_seq = zeros(1, N);       % Serial bit sequence (1..1023)
states_history = zeros(N, m); % History of all register states
for k = 1:N
    % 1. Save current state and output bit (cell 10)
    states_history(k, :) = state;
    prbs_seq(k) = state(10);
    % 2. Compute feedback bit (XOR between cell 10 and cell 7)
    feedback = xor(state(10), state(7));
    % 3. Right shift: feedback enters Q1, Q10 exits
    state = [feedback, state(1:end-1)];
end
% Verify results
fprintf('First 20 generated bits: ');
fprintf('%d', prbs_seq(1:20));
N0 = sum(prbs_seq == 0);
N1 = sum(prbs_seq == 1);
fprintf('\nNumber of 1s generated: %d (expected: %d)\n', N1, 2^(m-1));
fprintf('Number of 0s generated: %d (expected: %d)\n', N0, 2^(m-1) - 1);
% Verify that the state after 1023 cycles returns exactly to the initial seed
if isequal(state, states_history(1, :))
    disp('Maximum period verification: PASSED (returned to initial seed after 1023 cycles).');
end
% Conversion to bipolar sequence (0 -> -1, 1 -> +1) and plot
bipolar_seq = 2 * prbs_seq - 1;   
figure;
theme(gcf, 'light');
hold on;
stairs(1:20, bipolar_seq(1:20), 'LineWidth', 1); % square wave
grid on;
ylim([-1.5 1.5]);
xlim([0 21]);
xticks(1:20);
yticks([-1 0 1]);
xlabel('Bit index (n)');
ylabel('Bipolar Amplitude');
title('First 20 bits of the bipolar LFSR sequence (P(x) = x^{10} + x^7 + 1)');
%% --- Number of zero/ones, number of transitions ---
% make a table with N1 and N0
% Circular shift by 1 position to evaluate periodic transitions
shifted_seq = circshift(prbs_seq, -1);
transitions = xor(prbs_seq, shifted_seq);
NT  = sum(transitions == 1);
NNT = sum(transitions == 0);
T = table(N, N1, N0, NT, NNT, ...
    'VariableNames', {'Period_N', 'N1', 'N0', 'NT', 'NNT'});
disp(T);
% Verify ideal binary sequence balance property
% For an m-sequence (N = 2^m - 1), N1 and N0 differ by exactly 1 bit,
% so N0 ≈ N1 ≈ N/2 for large N.
if abs(N1 - N0) <= 1 && abs(NT - NNT) <= 1
    fprintf('Verified: N0 (%d) \x2248 N1 (%d) \x2248 N/2 (%.1f)\n', N0, N1, N/2);
    fprintf('Verified: NT (%d) \x2248 NNT (%d) \x2248 N/2 (%.1f)\n', NT, NNT, N/2);
else
    fprintf('Balance property NOT satisfied.\n');
end
%% --- Runs---
% 1. Align the sequence so that it does not start and end with the same bit
seq_runs = prbs_seq();      % to not lose the original value when rotate with circshift
shift_amount = 0;
while seq_runs(1) == seq_runs(end)
    seq_runs = circshift(seq_runs, 1);
    shift_amount = shift_amount + 1;
end
% 2. Count the length of each run using a simple loop
NR0 = zeros(1, m);   % zeros counter
NR1 = zeros(1, m);   % ones counter
run_lengths = [];
run_starts = [];
run_bits = [];
count = 1;
current_start = 1;
for k = 2:N
    if seq_runs(k) == seq_runs(k-1)
        count = count + 1;          % Bit equals previous bit: extend current run
    else
        run_lengths(end+1) = count;
        run_starts(end+1) = current_start;
        run_bits(end+1) = seq_runs(k-1);
        if seq_runs(k-1) == 0
            NR0(count) = NR0(count) + 1;
        else
            NR1(count) = NR1(count) + 1;
        end
        % Reset for the new run
        count = 1;
        current_start = k;
    end
end
% Save the last run of the sequence
run_lengths(end+1) = count;
run_starts(end+1)  = current_start;
run_bits(end+1)    = seq_runs(N);
if seq_runs(N) == 0
    NR0(count) = NR0(count) + 1;
else
    NR1(count) = NR1(count) + 1;
end
% 3. Build and print the summary table
NR_i = NR0 + NR1;
NR_T = sum(NR_i);
ideal_vals = (NR_T ./ (2 .^ (1:m)))';
fprintf('Total number of runs NR(T) = %d\n\n', NR_T);
T_runs = table((1:m)', NR0', NR1', NR_i', ideal_vals, ...
    'VariableNames', {'Length_i', 'NR0_i', 'NR1_i', 'NR_i', 'Ideal_NR_over_2i'});
disp(T_runs);
% 4. Plot: Run length vs. Starting point (different colors for 0 and 1)
idx_0 = (run_bits == 0);
idx_1 = (run_bits == 1);
figure;
theme(gcf, 'light');
stem(run_starts(idx_0), run_lengths(idx_0), 'Marker', 'none', 'DisplayName', 'Runs of 0s');
hold on;
stem(run_starts(idx_1), run_lengths(idx_1), 'Marker', 'none', 'DisplayName', 'Runs of 1s');
grid on;
xlim([1, N]);
ylim([0, m + 1]);
yticks(1:m);
xlabel('Starting index in the sequence (Starting point)');
ylabel('Run length (i)');
title('Run length vs. Starting point in the sequence (m = 10)');
legend('Location', 'best');

%% --- Periodic autocorrelation ---
% mean value E[ai] = 0
mean_val = mean(bipolar_seq);
fprintf('Mean of bipolar symbols E[a_i]: %.5f (\x2248 0)\n', mean_val);
% Range of tau required by the assignment: -(N-1)/2 <= tau <= +(N-1)/2
max_tau = (N - 1) / 2;            % Equals 511 for N = 1023
lags    = -max_tau : max_tau;     % Lag vector from -511 to +511
R_tau   = zeros(size(lags));
% For loop iterating over tau
for k = 1:length(lags)
    tau = lags(k);
    shifted_seq = circshift(bipolar_seq, tau);      % Compute periodic v'(i - tau)
    R_tau(k)    = sum(bipolar_seq .* shifted_seq);  % Sum from i=0 to N-1
end
% Plot autocorrelation
figure;
theme(gcf, 'light');
plot(lags, R_tau, 'LineWidth', 1.5);
grid on;
xlim([-max_tau, max_tau]);
ylim([-100, N + 100]);
xlabel('\tau');
ylabel('R_{v''}(\tau)');
title('Periodic Autocorrelation of the bipolar sequence (m = 10, N = 1023)');
%% --- truncated sequences MPSL ---
z_vals = 1:200;
for z = z_vals
    % Work directly on the bipolar sequence (-1, +1):
    bipolar_trunc = bipolar_seq(1 : end - z);
    N_trunc = length(bipolar_trunc); % Equals N - z
    % --- Compute properties on seq_trunc / bipolar_trunc ---
    % periodic autocorrelation on the truncated sequence
    max_tau = floor((N_trunc - 1) / 2);
    lags_trunc = -max_tau : max_tau;
    R_trunc = zeros(size(lags_trunc));
    for k = 1:length(lags_trunc)
        tau = lags_trunc(k);
        R_trunc(k) = sum(bipolar_trunc .* circshift(bipolar_trunc, tau));
    end
    % Metric to save for each z: maximum absolute value outside the peak
    MPSL(z) = max(abs(R_trunc(lags_trunc ~= 0)));
end
% PLOT 
figure;
theme(gcf, 'light');
plot(z_vals, MPSL, 'LineWidth', 1.5);
grid on;
ylim([min(MPSL) - 10, max(MPSL) + 10]);
xlabel('Number of removed bits (z)');
ylabel('Maximum Peak Side Lobe max_{\tau \neq 0} |R(\tau)|');
title('Autocorrelation degradation as a function of removed bits z');

%% --- truncated sequence with z = 23 ---
z = 23;
bipolar_seed_t = bipolar_seq(1 : N - z);
N_trunc = length(bipolar_seed_t); % Equals N - z

% Compute periodic autocorrelation on the truncated sequence
max_tau = floor((N_trunc - 1) / 2);
lags_trunc = -max_tau : max_tau;
R_trunc = zeros(size(lags_trunc));

for k = 1:length(lags_trunc)
    tau = lags_trunc(k);
    R_trunc(k) = sum(bipolar_seed_t .* circshift(bipolar_seed_t, tau));
end

% Maximum Peak Side Lobe (MPSL) for z = 23
MPSL_23 = max(abs(R_trunc(lags_trunc ~= 0)));
fprintf('For z = %d (length N_trunc = %d), MPSL = %d\n', z, N_trunc, MPSL_23);

% Plot periodic autocorrelation for z = 23
figure;
theme(gcf, "light");
plot(lags_trunc, R_trunc, 'LineWidth', 1.2);
grid on;
xlim([-max_tau, max_tau]);
ylim([-100, N_trunc + 100]);
xlabel('\tau');
ylabel('R_{v''_{23}}(\tau)');
seed_str = sprintf('%d', state);
title({sprintf('Periodic Autocorrelation of Truncated Sequence (z = %d, MPSL = %d)', z, MPSL_23), ...
    sprintf('Seed = [%s]', seed_str)});

% repeat and plot with another starting seed
seed_23 = [1 1 0 0 1 1 0 0 1 1];
prbs_seq_23 = zeros(1, N);
seed_history = zeros(N, m); % History of all register states

for k = 1:N
    % 1. Save current state and output bit (cell 10)
    seed_history(k, :) = seed_23;
    prbs_seq_23(k) = seed_23(10);
    % 2. Compute feedback bit (XOR between cell 10 and cell 7)
    feedback = xor(seed_23(10), seed_23(7));
    % 3. Right shift: feedback enters Q1, Q10 exits
    seed_23 = [feedback, seed_23(1:end-1)];
end

% Conversion to bipolar sequence for the new seed
bipolar_seq_23 = 2 * prbs_seq_23 - 1;

% Verify results
fprintf('First 20 generated bits: ');
fprintf('%d', prbs_seq_23(1:20));
N0 = sum(prbs_seq_23 == 0);
N1 = sum(prbs_seq_23 == 1);
fprintf('\nNumber of 1s generated: %d (expected: %d)\n', N1, 2^(m-1));
fprintf('Number of 0s generated: %d (expected: %d)\n', N0, 2^(m-1) - 1);
% Verify that the state after 1023 cycles returns exactly to the initial seed
if isequal(seed_23, seed_history(1, :))
    disp('Maximum period verification: PASSED (returned to initial seed after 1023 cycles).');
end

% z = 23
bipolar_seed_t = bipolar_seq_23(1 : N - z);
N_trunc = length(bipolar_seed_t); % Equals N - z

% Compute periodic autocorrelation on the truncated sequence
max_tau = floor((N_trunc - 1) / 2);
lags_trunc = -max_tau : max_tau;
R_trunc = zeros(size(lags_trunc));

for k = 1:length(lags_trunc)
    tau = lags_trunc(k);
    R_trunc(k) = sum(bipolar_seed_t .* circshift(bipolar_seed_t, tau));
end

% Maximum Peak Side Lobe (MPSL) for z = 23
MPSL_23 = max(abs(R_trunc(lags_trunc ~= 0)));
fprintf('For z = %d (length N_trunc = %d), MPSL = %d\n', z, N_trunc, MPSL_23);

% Plot periodic autocorrelation for z = 23
figure;
theme(gcf, "light");
plot(lags_trunc, R_trunc, 'LineWidth', 1.2);
grid on;
xlim([-max_tau, max_tau]);
ylim([-100, N_trunc + 100]);
xlabel('\tau');
ylabel('R_{v''_{23}}(\tau)');
seed_str = sprintf('%d', seed_23); % seed [1 1 0 0 1 1 0 0 1 1] in '1100110011'
title({sprintf('Periodic Autocorrelation of Truncated Sequence (z = %d, MPSL = %d)', z, MPSL_23), ...
    sprintf('Seed = [%s]', seed_str)});