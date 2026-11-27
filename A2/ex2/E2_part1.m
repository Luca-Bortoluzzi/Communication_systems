clc;
clear;
close all;

%% -- Exercise 2 - CDMA ---
%% --- Part 1: 2-PAM ---
% Plot the analytical BER vs Eb/N0 from 0 to 12 db

Eb = 1;     % energy per bit (pulse unitary energy)

EbN0_db = 0:12;     % Eb / N0
EbN0_lin = 10.^(EbN0_db / 10);  % Linear Conversion
P_a = 1/2 * erfc(sqrt(EbN0_lin));

figure;
theme(gcf, "light");
semilogy(EbN0_db, P_a, 'LineWidth', 2);
xlabel("Eb / N0 (dB unit)");
ylabel("BER");
title("Analytical BER of 2-PAM over a AWGN channel");
axis("tight");
grid on;
exportgraphics(gcf, 'A2_2_1.png', 'Resolution', 300);

%% --- Monte Carlo simulation (0 to 8 dB) ---
EbN0_sim_db = 0:8;
EbN0_sim_lin = 10.^(EbN0_sim_db / 10);
BER_sim = zeros(size(EbN0_sim_db));

for i = 1:length(EbN0_sim_db)
    errors = 0;
    bits = 0;

    N0 = Eb / EbN0_sim_lin(i);
    sigma = sqrt(N0 / 2);

    % Repeat until 100 errors
    while errors < 100
        b = randi([0, 1]);
        bipolar = 2*b - 1;          % 0 -> -1, 1 -> +1

        noise = sigma * randn;      % AWGN ~ N(0, N0/2)
        r = bipolar + noise;

        % Minimum distance detector (0 threshold)
        if r < 0
            a = 0;
        else
            a = 1;
        end

        if a ~= b
            errors = errors + 1;
        end
        bits = bits + 1;
    end

    BER_sim(i) = errors / bits;
end

%% --- Plot Analytic vs Simulated ---
figure;
theme(gcf, "light");
semilogy(EbN0_db, P_a, 'b-', 'LineWidth', 2); hold on;
semilogy(EbN0_sim_db, BER_sim, 'ro--', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel("E_b / N_0 (dB)");
ylabel("BER");
title("2-PAM: Analytical vs Simulated BER");
legend("Analytical", "Monte Carlo Simulation");
grid on;
exportgraphics(gcf, 'A2_2_2.png', 'Resolution', 300);

%% --- comment on the results

% The monte Carlo simulation on 100 errors (red curve) follows almost perfectly the
% Analytical curve (blue curve) with a little misalignment on 2dB
