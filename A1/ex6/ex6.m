clc;
close all;
clear;

%% --- 1. LFSR Randomizer Design and m-sequence Generation ---
% Number of shift-register cells (flip-flops)
m = 10;
% Primitive generator polynomial: P(x) = x^{10} + x^7 + 1

% Initial non-zero state (Seed 1)
seed_1 = [1 1 1 1 1 0 0 0 0 0];

% Generate one complete period of the m-sequence (N = 2^m - 1 = 1023 bits)
N = 2^m - 1;
prbs_s1 = zeros(1, N);        % Output serial binary sequence (0s and 1s)
seed_1_history = zeros(N, m); % History of all internal register states

for k = 1:N
    % 1. Store current register state and output bit from cell 10 (Q10)
    seed_1_history(k, :) = seed_1;
    prbs_s1(k) = seed_1(10);
    
    % 2. Compute feedback bit via XOR between taps Q10 and Q7 (x^10 + x^7 + 1)
    feedback = xor(seed_1(10), seed_1(7));
    
    % 3. Shift register right: feedback enters Q1, Q10 exits
    seed_1 = [feedback, seed_1(1:end-1)];
end

% Verify generated bit counts and maximum-period property
fprintf('First 20 generated bits: ');
fprintf('%d', prbs_s1(1:20));
N0 = sum(prbs_s1 == 0);
N1 = sum(prbs_s1 == 1);
fprintf('\nNumber of 1s generated: %d (expected: %d)\n', N1, 2^(m-1));
fprintf('Number of 0s generated: %d (expected: %d)\n', N0, 2^(m-1) - 1);

% Verify that the register returns to the initial seed after exactly N = 1023 cycles
if isequal(seed_1, seed_1_history(1, :))
    disp('Maximum period verification: PASSED (returned to initial seed after 1023 cycles).');
end

% Map binary symbols to bipolar levels (0 -> -1, 1 -> +1) and plot first 20 bits
bipolar_s1 = 2 * prbs_s1 - 1;   
figure;
theme(gcf, 'light');
hold on;
stairs(1:20, bipolar_s1(1:20), 'LineWidth', 1); % Square-wave representation
grid on;
ylim([-1.5 1.5]);
xlim([0 21]);
xticks(1:20);
yticks([-1 0 1]);
xlabel('Bit index (n)');
ylabel('Bipolar Amplitude');
title('First 20 bits of the bipolar LFSR sequence (P(x) = x^{10} + x^7 + 1)');

%% --- 2. Balance and Transition Properties ---
% Circularly shift by 1 position to evaluate periodic bit transitions
shifted_seq = circshift(prbs_s1, -1);
transitions = xor(prbs_s1, shifted_seq);
NT  = sum(transitions == 1); % Number of transitions (0->1 or 1->0)
NNT = sum(transitions == 0); % Number of non-transitions (0->0 or 1->1)

% Display summary table of bit and transition counts
T = table(N, N1, N0, NT, NNT, ...
    'VariableNames', {'Period_N', 'N1', 'N0', 'NT', 'NNT'});
disp(T);

% Verify ideal random binary sequence balance properties:
% For an m-sequence (N = 2^m - 1), N1 and N0 differ by 1 bit,
% so N0 ≈ N1 ≈ N/2 and NT ≈ NNT ≈ N/2 for large N.
if abs(N1 - N0) <= 1 && abs(NT - NNT) <= 1
    fprintf('Verified: N0 (%d) \x2248 N1 (%d) \x2248 N/2 (%.1f)\n', N0, N1, N/2);
    fprintf('Verified: NT (%d) \x2248 NNT (%d) \x2248 N/2 (%.1f)\n', NT, NNT, N/2);
else
    fprintf('Balance property NOT satisfied.\n');
end

%% --- 3. Run-Length Distribution Analysis ---
% 1. Circularly align the sequence so that the first and last bits differ,
%    preventing boundary-wrapping runs from being split into two shorter runs.
seq_runs = prbs_s1();
shift_amount = 0;
while seq_runs(1) == seq_runs(end)
    seq_runs = circshift(seq_runs, 1);
    shift_amount = shift_amount + 1;
end

% 2. Measure the length and starting index of every run of 0s and 1s
NR0 = zeros(1, m); % Run-length histogram for 0s (lengths 1..m)
NR1 = zeros(1, m); % Run-length histogram for 1s (lengths 1..m)
run_lengths = [];
run_starts = [];
run_bits = [];
count = 1;
current_start = 1;

for k = 2:N
    if seq_runs(k) == seq_runs(k-1)
        count = count + 1; % Current bit matches previous bit: extend run
    else
        % Save completed run statistics
        run_lengths(end+1) = count;
        run_starts(end+1) = current_start;
        run_bits(end+1) = seq_runs(k-1);
        if seq_runs(k-1) == 0
            NR0(count) = NR0(count) + 1;
        else
            NR1(count) = NR1(count) + 1;
        end
        % Reset counters for the new run
        count = 1;
        current_start = k;
    end
end

% Record the final run at the end of the aligned period
run_lengths(end+1) = count;
run_starts(end+1)  = current_start;
run_bits(end+1)    = seq_runs(N);
if seq_runs(N) == 0
    NR0(count) = NR0(count) + 1;
else
    NR1(count) = NR1(count) + 1;
end

% 3. Build and display summary table comparing measured runs with ideal values (NR_T / 2^i)
NR_i = NR0 + NR1;
NR_T = sum(NR_i);
ideal_vals = (NR_T ./ (2 .^ (1:m)))';
fprintf('Total number of runs NR(T) = %d\n\n', NR_T);
T_runs = table((1:m)', NR0', NR1', NR_i', ideal_vals, ...
    'VariableNames', {'Length_i', 'NR0_i', 'NR1_i', 'NR_i', 'Ideal_NR_over_2i'});
disp(T_runs);

% 4. Plot run length vs. starting index (separate series for 0-runs and 1-runs)
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

%% --- 4. Periodic Autocorrelation of the Full m-sequence ---
% Compute sample mean of bipolar symbols (expected E[a_i] ≈ 0)
mean_val = mean(bipolar_s1);
fprintf('Mean of bipolar symbols E[a_i]: %.5f (\x2248 0)\n', mean_val);

% Evaluate periodic autocorrelation R_{v'}(tau) over -(N-1)/2 <= tau <= +(N-1)/2
max_tau = (N - 1) / 2;            % Equals 511 for N = 1023
lags    = -max_tau : max_tau;     % Lag vector from -511 to +511
R_tau   = zeros(size(lags));

for k = 1:length(lags)
    tau = lags(k);
    shifted_seq = circshift(bipolar_s1, tau);      % Circularly shifted sequence v'(i - tau)
    R_tau(k)    = sum(bipolar_s1 .* shifted_seq);  % Periodic correlation sum
end

% Plot periodic autocorrelation (two-level: N at tau = 0, -1 for tau != 0)
figure;
theme(gcf, 'light');
plot(lags, R_tau, 'LineWidth', 1.5);
grid on;
xlim([-max_tau, max_tau]);
ylim([-100, N + 100]);
xlabel('\tau');
ylabel('R_{v''}(\tau)');
title('Periodic Autocorrelation of the bipolar sequence (m = 10, N = 1023)');

%% --- 5. Truncated Sequences: MPSL Degradation Analysis (1 <= z <= 200) ---
z_vals = 1:200;
for z = z_vals
    % Truncate the last z symbols from the bipolar sequence (new length = N - z)
    bipolar_trunc = bipolar_s1(1 : end - z);
    N_trunc = length(bipolar_trunc);
    
    % Compute periodic autocorrelation on the truncated sequence
    max_tau = floor((N_trunc - 1) / 2);
    lags_trunc = -max_tau : max_tau;
    R_trunc = zeros(size(lags_trunc));
    for k = 1:length(lags_trunc)
        tau = lags_trunc(k);
        R_trunc(k) = sum(bipolar_trunc .* circshift(bipolar_trunc, tau));
    end
    
    % Extract Maximum Peak Side Lobe: MPSL(z) = max_{tau != 0} |R(tau)|
    MPSL(z) = max(abs(R_trunc(lags_trunc ~= 0)));
end

% Plot MPSL(z) versus number of removed bits z
figure;
theme(gcf, 'light');
plot(z_vals, MPSL, 'LineWidth', 1.5);
grid on;
ylim([min(MPSL) - 10, max(MPSL) + 10]);
xlabel('Number of removed bits (z)');
ylabel('Maximum Peak Side Lobe max_{\tau \neq 0} |R(\tau)|');
title('Autocorrelation degradation as a function of removed bits z');

%% --- 6. Truncated Sequence with z = 23: Comparison of Two Starting Seeds ---
% --- Case A: Seed 1 = [1 1 1 1 1 0 0 0 0 0] ---
z = 23;
bipolar_trunc_s1 = bipolar_s1(1 : N - z); % Truncate last 23 symbols (length = 1000)
N_trunc = length(bipolar_trunc_s1);

% Compute periodic autocorrelation for the truncated sequence (Seed 1)
max_tau = floor((N_trunc - 1) / 2);
lags_trunc = -max_tau : max_tau;
R_trunc_s1 = zeros(size(lags_trunc));

for k = 1:length(lags_trunc)
    tau = lags_trunc(k);
    R_trunc_s1(k) = sum(bipolar_trunc_s1 .* circshift(bipolar_trunc_s1, tau));
end

% Compute MPSL for z = 23 (Seed 1)
MPSL_1= max(abs(R_trunc_s1(lags_trunc ~= 0)));
fprintf('For z = %d (length N_trunc = %d), MPSL = %d\n', z, N_trunc, MPSL_1);

% Plot periodic autocorrelation of the truncated sequence (Seed 1)
figure;
theme(gcf, "light");
plot(lags_trunc, R_trunc_s1, 'LineWidth', 1.2);
grid on;
xlim([-max_tau, max_tau]);
ylim([-100, N_trunc + 100]);
xlabel('\tau');
ylabel('R_{v''_{23}}(\tau)');
seed_str = sprintf('%d', seed_1);
title({sprintf('Periodic Autocorrelation of Truncated Sequence (z = %d, MPSL = %d)', z, MPSL_1), ...
    sprintf('Seed = [%s]', seed_str)});

% --- Case B: Seed 2 = [1 1 0 0 1 1 0 0 1 1] ---
seed_2 = [1 1 0 0 1 1 0 0 1 1];
prbs_s2 = zeros(1, N);
seed_history = zeros(N, m); % History of all internal register states for Seed 2

% Generate complete m-sequence starting from Seed 2
for k = 1:N
    % 1. Store current register state and output bit from cell 10 (Q10)
    seed_history(k, :) = seed_2;
    prbs_s2(k) = seed_2(10);
    
    % 2. Compute feedback bit via XOR between taps Q10 and Q7 (x^10 + x^7 + 1)
    feedback = xor(seed_2(10), seed_2(7));
    
    % 3. Shift register right: feedback enters Q1, Q10 exits
    seed_2 = [feedback, seed_2(1:end-1)];
end

% Convert Seed 2 binary sequence to bipolar format (0 -> -1, 1 -> +1)
bipolar_s2 = 2 * prbs_s2 - 1;

% Truncate last z = 23 symbols (length = 1000)
bipolar_trunc_s2 = bipolar_s2(1 : N - z);
N_trunc = length(bipolar_trunc_s2);

% Compute periodic autocorrelation for the truncated sequence (Seed 2)
max_tau = floor((N_trunc - 1) / 2);
lags_trunc = -max_tau : max_tau;
R_trunc_s2 = zeros(size(lags_trunc));

for k = 1:length(lags_trunc)
    tau = lags_trunc(k);
    R_trunc_s2(k) = sum(bipolar_trunc_s2 .* circshift(bipolar_trunc_s2, tau));
end

% Compute MPSL for z = 23 (Seed 2)
MPSL_2 = max(abs(R_trunc_s2(lags_trunc ~= 0)));
fprintf('For z = %d (length N_trunc = %d), MPSL = %d\n', z, N_trunc, MPSL_2);

% Plot periodic autocorrelation of the truncated sequence (Seed 2)
figure;
theme(gcf, "light");
plot(lags_trunc, R_trunc_s2, 'LineWidth', 1.2);
grid on;
xlim([-max_tau, max_tau]);
ylim([-100, N_trunc + 100]);
xlabel('\tau');
ylabel('R_{v''_{23}}(\tau)');
seed_str = sprintf('%d', seed_2);
title({sprintf('Periodic Autocorrelation of Truncated Sequence (z = %d, MPSL = %d)', z, MPSL_2), ...
    sprintf('Seed = [%s]', seed_str)});