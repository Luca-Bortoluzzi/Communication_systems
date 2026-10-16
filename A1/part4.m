clc;
clear;
close all;

%% part 4
% previous part to make work the code
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

% 4.6 (ppt part)
figure;
semilogy(p, PUE(1, :), '-o', 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'Target: 10 Errori');
hold on;
semilogy(p, PUE(2, :), '-s', 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'Target: 100 Errori');
hold on;
% Theoretical Analytic Curve given by p_interval values
semilogy(p, PUE_analytical, '--r', 'LineWidth', 2.0, 'DisplayName', 'Theoretical Analytic Curve');

grid on;
% to see how much the curves are close to the theoretical remove 'ylim'
ylim([min(PUE(:))*0.5, max(PUE(:))*2]);     % add limits and zoom the y axis to better reading

xlabel('Error Probability in the channel (p)');
ylabel('Probability of Undetected Error (PUE)');
title('Comparison between Monte Carlo simulations and Analytic Curve');
legend('Location', 'best');
set(gca, 'XDir', 'reverse'); % Invert the X axis to better see p from 0.10 to 0.06