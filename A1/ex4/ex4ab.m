% Exercise 4 A1

clc;
close all;
clear;

%% 4a 
% Verify the code is perfect

G = [1 0 0 0 1 1 0;
    0 1 0 0 1 0 1;
    0 0 1 0 0 1 1;
    0 0 0 1 1 1 1];

% C(7, 4)
n = size(G, 2);     % taking the number of columns
k = size(G, 1);     % taking the number of rows

v = dec2bin(0:2^k-1) - '0';  

% Matrix multiplication modulo 2 generates all 2^k possible transmitted Codewords.
c = mod(v*G, 2);   

% Let's calculate the hamming weight with zeros
Hw_0 = sum(c, 2);

% and without zeros
Hw = Hw_0(Hw_0 > 0);

% minimum distance excluding zero values
dmin = min(Hw);

% we also define the error correction capability as
t = (dmin - 1) / 2;

% A code is perfect if the Hamming spheres of radius t exactly cover the received-vector space without overlap.
% D(c) = S(c, t) for c in C         [perfect codes]

% we define S(c, t)
S = 0;
for i=0:t
    S = S + nchoosek(n, i);
end

% condition for perfect codes
condition = 2^k * S;

if(condition == 2^n)
    fprintf("This is a perfect code.\n")
else
    fprintf("This isn't a perfect code.\n");
end

% Assuming transmission over a BSC with crossover probability p and minimum Hamming distance
% decoding, derive the analytical expression for the word error probability PW(E) as a function of p.

p = logspace(-4, log10(0.5), 100);

p_corr = zeros(size(p));

for i=0:t
    p_corr = p_corr + nchoosek(n, i) .* p.^i .* (1-p).^(n - i);
end

p_word_err = 1 - p_corr;

figure;
theme(gcf, "light"); % Imposta l'intera figura (sfondo, assi, scritte, legenda) come nel tema chiaro
loglog(p, p_word_err, 'b-', 'LineWidth', 2);
grid on;
xlabel('Bit error probability in channel (p)');
ylabel('Codeword error probability P_W(E)');
title(sprintf('Theoretical P_W(E) for a (%d, %d) perfect code (t = %d)', n, k, t));
legend('Theoretical P_W(E)', 'Location', 'southeast');
axis('tight');

%% 4b

% Ex1 part 4 Monte Carlo simulation - Word Error Probability PW(E)
target_errors = [10, 100];     % Stop simulation upon reaching 10 or 100 word errors
rng(70);    % fixed seed
p = [0.10, 0.09, 0.08, 0.07, 0.06];
PW_E = zeros(length(target_errors), length(p));
PW_E_analytical = zeros(1, length(p));

% --- 1. MONTE CARLO SIMULATION ---
for j = 1:length(target_errors)
    for i = 1:length(p)
        word_errors = 0;        % Reset error counter for each value of p
        total_transmitted = 0;

        while word_errors < target_errors(j)
            % 4.1 Message generation and transmission over BSC
            % the codeword generated and random message doesn't influence
            % the probability word error, so we don't need them

            e = rand(1, n) < p(i);              % BSC noise pattern generation

            % 4.2 Word error condition under Minimum Distance Decoding
            % A word error occurs whenever the number of channel errors exceeds t
            num_errors = sum(e);
            if num_errors > t
                word_errors = word_errors + 1;
            end

            total_transmitted = total_transmitted + 1;
        end

        % Empirical Word Error Probability
        PW_E(j, i) = word_errors / total_transmitted;
    end
end

% --- 2. ANALYTICAL CALCULATION ---
for i = 1:length(p)
    p_val = p(i);
    pwe_sum = 0;
    % Binomial tail summation from w = t + 1 to n (more than t errors)
    %        All error patterns of weight at most t are corrected.
    %        No error pattern of weight greater than t is corrected.
    % so we use t+1
    for w = (t + 1):n       
        pwe_sum = pwe_sum + nchoosek(n, w) * (p_val^w) * ((1 - p_val)^(n - w));
    end
    PW_E_analytical(i) = pwe_sum;
end

% --- 3. COMPARISON PLOT ---
figure;
theme(gcf, "light");
loglog(p, PW_E_analytical, 'b-', 'LineWidth', 2, 'DisplayName', 'Analytical'); hold on;
loglog(p, PW_E(1, :), 'r--o', 'LineWidth', 1.5, 'DisplayName', sprintf('MC (target = %d)', target_errors(1)));
loglog(p, PW_E(2, :), '--s', 'LineWidth', 1.5, 'Color', '#77AC30', 'DisplayName', sprintf('MC (target = %d)', target_errors(2)));
grid on;
xlabel('Channel bit error probability (p)');
ylabel('Word error probability P_W(E)');
title('Monte Carlo vs Analytical Comparison for P_W(E)');
legend('Location', 'southeast');