clc;
clear;
close all;

%% part 4
% previous part to make the code work
A = [1 1 1 0 0;
    0 1 1 1 0;
    0 0 1 1 1;
    1 0 0 1 1;
    1 1 0 0 1;];

k = 5;   % information length != codeword length = (info length + redundancy bits)
G = [eye(5) A];
v = dec2bin(0:2^k-1) - '0';
c1 = mod(v*G, 2);
hw=zeros(32, 1);
mult=zeros(11,1);

for i=1:32
    hw(i, 1) = sum(c1(i, : ));
    a1 = hw(i, 1);
    mult(a1+1, 1) = mult(a1+1, 1) + 1;
end

% real start of the part4

p = [0.10, 0.09, 0.08, 0.07, 0.06];     % given prob for part 4

n = 10; % codeword length
PUE_analytical = zeros(1, length(p));

% BSC modeling: y = c_rand + e | where y is the codeword received and "e" are errors during the transmission

P = G(:, k+1:end);     %parity submatrix (from 6 to 10 columns [5x5])
H = [P', eye(n-k)];    % Syndrome matrix

target_errors = [10, 100];     % stop the cycle at undetected error = 10;
totalError = 0;
PUE = zeros(length(target_errors), length(p));
% we need to stop the cycle at 10 undetected errors and then at 100
for j=1:2
    for i=1:length(p)

        undetected_errors = 0;  % self-explanatory
        total_transmitted = 0;  % //

        while undetected_errors < target_errors(j)
            % 4.1
            u = randi([0, 1], 1, k);     % random generation of a [1 x k] only for 0 and 1
            c_rand = mod(u * G, 2);    % random gen of a vector 1 x 10 (codeword length)

            e = rand(1, length(c_rand)) < p(i);  % generate random noise in BSC with less probability than current
            y = mod(c_rand + e, 2);              % received message

            % 4.2 Syndrome calculus
            s = mod(y * H', 2);

            % 4.3 Error detection
            if all(s == 0) && ~isequal(y, c_rand)
                undetected_errors = undetected_errors + 1;  % undetected increment
            end

            total_transmitted = total_transmitted + 1;  % transmitted anyway

        end
        PUE(j, i) = undetected_errors / total_transmitted;
    end
end

for i=1:length(p)
    pue_sum = 0;
    p_val = p(i);

    for w=1:n
        pue_sum = pue_sum + mult(w+1) * (p_val^w) * ((1 - p_val)^(n - w));
    end
    PUE_analytical(i) = pue_sum;
end

% Colori coordinati con la presentazione Beamer
plot_blue   = [30, 136, 229] / 255;
plot_orange = [251, 140, 0] / 255;
plot_red    = [229, 57, 53] / 255;

fig_p4 = figure('Color', 'w', 'Position', [100, 100, 760, 400]);

semilogy(p, PUE(1, :), '-o', 'Color', plot_blue, 'MarkerFaceColor', plot_blue, ...
    'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', 'Target: 10 Errors');
hold on;
semilogy(p, PUE(2, :), '-s', 'Color', plot_orange, 'MarkerFaceColor', plot_orange, ...
    'LineWidth', 1.8, 'MarkerSize', 6, 'DisplayName', 'Target: 100 Errors');
semilogy(p, PUE_analytical, '--', 'Color', plot_red, ...
    'LineWidth', 2.2, 'DisplayName', 'Theoretical Analytic Curve');

grid on; box off;
ylim([min([PUE(:); PUE_analytical(:)])*0.55, max([PUE(:); PUE_analytical(:)])*1.8]);
xlim([0.058, 0.102]);
xticks(fliplr(p)); % Mantiene i tick ordinati: 0.06, 0.07, 0.08, 0.09, 0.10
xlabel('BSC crossover probability p', 'FontSize', 11);
ylabel('Probability of Undetected Error (P_{UE})', 'FontSize', 11);
title('Monte Carlo Simulation vs. Analytical P_{UE}', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northwest', 'FontSize', 10);
exportgraphics(fig_p4, 'E1_P4_Q1.png', 'Resolution', 300);
