%% Exercise 2 - Part 2: Single-user DSSS
clc;
clear;
close all;

%% --- System Parameters ---
N = 12;                         % Spreading factor (number of chips per bit)
Eb = N;                         % Energy per bit (N chips of unit energy)
target_errors = 100;            % Minimum number of bit errors to collect per SNR point

EbN0_db = 0:8;                  % Eb/N0 range in dB (step of 1 dB)
EbN0_lin = 10.^(EbN0_db / 10);  % Convert Eb/N0 from dB to linear scale

% Analytical BER curve (equivalent to 2-PAM on AWGN channel)
P_a = 1/2 * erfc(sqrt(EbN0_lin));

%% --- Spreading Sequence Generation (User U1) ---
rng(42);                        % Set random seed for reproducibility
c1 = randi([0, 1], 1, N);       % Generate a fixed random binary sequence of length N
c1_bipolar = 2*c1 - 1;          % Bipolar mapping: 0 -> -1, 1 -> +1

%% --- Monte Carlo Simulation ---
BER_sim = zeros(size(EbN0_db)); % Preallocate array for simulated BER values

for i = 1:length(EbN0_db)
    errors = 0;                 % Reset error counter for current Eb/N0
    bits = 0;                   % Reset transmitted bit counter

    N0 = Eb / EbN0_lin(i);      % Compute noise power spectral density N0 (with Eb = N)
    sigma = sqrt(N0 / 2);       % Noise standard deviation per chip (variance = N0/2)

    % Run simulation until at least 100 bit errors are observed
    while errors < target_errors
        % 1. Generate a random information bit and map it to bipolar (-1 or +1)
        b = randi([0, 1]);
        bipolar = 2*b - 1;

        % 2. Spreading: multiply the bipolar symbol by the N-chip spreading sequence
        s = bipolar * c1_bipolar;

        % 3. AWGN Channel: add N independent Gaussian noise samples
        noise = sigma * randn(1, N);
        r = s + noise;

        % 4. Despreading: project received vector r onto c1_bipolar.
        % dot(A, B) computes the scalar (inner) product sum(A .* B), correlating
        % the received chips with the known user sequence to accumulate bit energy.
        y = dot(r, c1_bipolar);

        % 5. Minimum-distance decision (sign threshold at 0)
        if y < 0
            a = 0;
        else
            a = 1;
        end

        % Count bit errors and total simulated bits
        if a ~= b
            errors = errors + 1;
        end
        bits = bits + 1;
    end

    % Estimate BER for the current Eb/N0 value
    BER_sim(i) = errors / bits;
end

%% --- Plot Results ---
figure;
theme(gcf, "light");
semilogy(EbN0_db, P_a, 'b-', 'LineWidth', 2); hold on;
semilogy(EbN0_db, BER_sim, 'ro--', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel("E_b / N_0 (dB)");
ylabel("Bit Error Probability P_b(e)");
title("Single-user DSSS (N = 12): Analytical vs Simulated BER");
legend("Analytical (2-PAM)", "Simulated DSSS");
grid on;

%% --- comment on the results ---

% Comment on Single-user DSSS results (Point 6):
% The simulated BER curve for the single-user DSSS system closely matches the analytical 2-PAM curve $P_b(e) = \frac{1}{2}\operatorname{erfc}(\sqrt{E_b/N_0})$, 
% showing the same performance observed in Part 1. Since only one user is active and there is no multi-user interference, 
% spreading each bit over $N = 12$ unit-energy chips ($E_b = N$) and despreading at the receiver corresponds to 
% projecting the received signal onto the normalized basis vector $\underline{c}_1'/\sqrt{N}$. Consequently, at the same $E_b/N_0$, 
% spreading neither improves nor degrades the error rate over an AWGN channel; its role is rather to enable multi-user separation in CDMA. 
% Minor fluctuations at low $E_b/N_0$ ($0\text{–}3\text{ dB}$) are solely due to the statistical variance of the Monte Carlo estimation when stopping at $100$ bit errors.   