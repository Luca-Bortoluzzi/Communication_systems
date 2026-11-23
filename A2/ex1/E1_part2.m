clc; 
clear;
close all;

%% A2.1.5 (2 Gold seq with random i >= 0)
m = 7;          % number of cells (LFSR)
L = 2^m - 1;    % sequence length (127)
lags = 0 : (L-1); % symmetric lag axis for fftshift

rng(32); % 109 e 47 come valori random (TOGLIERE SUCCESSIVAMENTE)

% Valid indices for comm.GoldSequence with m=7 are [-2, L-2] -> [0, 125] for i >= 0
i = randi([0, L - 2], 2, 1);

goldseq = comm.GoldSequence('FirstPolynomial','x^7+x^4+1',...
    'SecondPolynomial','x^7+x+1',...
    'FirstInitialConditions',[0 0 0 0 0 0 1],...
    'SecondInitialConditions',[0 0 0 0 0 0 1],...
    'Index', i(1), 'SamplesPerFrame', L);

x = goldseq()';
xb = 1 - 2*x;

release(goldseq);          % Sblocca l'oggetto
goldseq.Index = i(2);      % Cambia solo l'indice
x1 = goldseq()';           % Genera la seconda sequenza
x1b = 1 - 2*x1;
%% A2.1.2 (Plot Cross-correlation)

R_cross = ifft(fft(xb).*conj(fft(x1b))); % non-normalized periodic autocorr.
R_cross = fftshift(R_cross);

figure;
theme(gcf, "light");
plot(lags, R_cross);
grid on;
xlabel('Lag');
ylabel('Cross-correlation');
title(sprintf("Cross-correlation of two Gold codes with N = %d ", L));
axis('tight');

%% A2.1.7 Generate two m-seq with diff polyms
N = L;

% --- First m-sequence: P1(x) = x^7 + x^4 + 1 ---
seed_1 = [1 1 1 1 0 0 0];
prbs_s1 = zeros(1, N);        % Output serial binary sequence (0s and 1s)
seed_1_history = zeros(N, m); % History of all internal register states

for k = 1:N
    % 1. Store current register state and output bit from cell 7 (Q7)
    seed_1_history(k, :) = seed_1;
    prbs_s1(k) = seed_1(7);

    % 2. Compute feedback bit via XOR between taps Q7 and Q4 (x^7 + x^4 + 1)
    feedback_1 = xor(seed_1(7), seed_1(4));

    % 3. Shift register right: feedback enters Q1, Q7 exits
    seed_1 = [feedback_1, seed_1(1:end-1)];
end

% Verify generated bit counts and maximum-period property for m-seq 1
fprintf('--- m-sequence 1 (x^7 + x^4 + 1) ---\n');
fprintf('First 20 generated bits: ');
fprintf('%d', prbs_s1(1:20));
N0_1 = sum(prbs_s1 == 0);
N1_1 = sum(prbs_s1 == 1);
fprintf('\nNumber of 1s generated: %d (expected: %d)\n', N1_1, 2^(m-1));
fprintf('Number of 0s generated: %d (expected: %d)\n', N0_1, 2^(m-1) - 1);

if isequal(seed_1, seed_1_history(1, :))
    fprintf('Maximum period verification: PASSED (returned to initial seed after %d cycles).\n\n', N);
end

% --- Second m-sequence: P2(x) = x^7 + x^2 + 1 ---
seed_2 = [1 1 1 1 0 0 0];
prbs_s2 = zeros(1, N);        % Output serial binary sequence (0s and 1s)
seed_2_history = zeros(N, m); % History of all internal register states

for k = 1:N
    % 1. Store current register state and output bit from cell 7 (Q7)
    seed_2_history(k, :) = seed_2;
    prbs_s2(k) = seed_2(7);

    % 2. Compute feedback bit via XOR between taps Q7 and Q1 (x^7 + x + 1)
    feedback_2 = xor(seed_2(7), seed_2(2));

    % 3. Shift register right: feedback enters Q1, Q7 exits
    seed_2 = [feedback_2, seed_2(1:end-1)];
end

% Verify generated bit counts and maximum-period property for m-seq 2
fprintf('--- m-sequence 2 (x^7 + x^2 + 1) ---\n');
fprintf('First 20 generated bits: ');
fprintf('%d', prbs_s2(1:20));
N0_2 = sum(prbs_s2 == 0);
N1_2 = sum(prbs_s2 == 1);
fprintf('\nNumber of 1s generated: %d (expected: %d)\n', N1_2, 2^(m-1));
fprintf('Number of 0s generated: %d (expected: %d)\n', N0_2, 2^(m-1) - 1);

if isequal(seed_2, seed_2_history(1, :))
    fprintf('Maximum period verification: PASSED (returned to initial seed after %d cycles).\n\n', N);
end

% Map binary symbols to bipolar levels (0 -> +1, 1 -> -1) consistently
bipolar_s1 = 1 - 2 * prbs_s1;
bipolar_s2 = 1 - 2 * prbs_s2;

% Plot first 20 bits of the first m-sequence
figure;
theme(gcf, 'light');
stairs(1:20, bipolar_s1(1:20), 'LineWidth', 2); % Square-wave representation
grid on;
ylim([-1.5 1.5]);
xlim([0 21]);
xticks(1:20);
yticks([-1 0 1]);
xlabel('Bit index (n)');
ylabel('Bipolar Amplitude');
title('First 20 bits of the bipolar LFSR sequence (P_1(x) = x^{7} + x^4 + 1)');

% Plot first 20 bits of the second m-sequence
figure;
theme(gcf, 'light');
stairs(1:20, bipolar_s2(1:20), 'LineWidth', 2); % Square-wave representation
grid on;
ylim([-1.5 1.5]);
xlim([0 21]);
xticks(1:20);
yticks([-1 0 1]);
xlabel('Bit index (n)');
ylabel('Bipolar Amplitude');
title('First 20 bits of the bipolar LFSR sequence (P_2(x) = x^{7} + x^2 + 1)');

%% A2.1.2.8

R_cross = ifft(fft(bipolar_s1).*conj(fft(bipolar_s2))); % non-normalized periodic autocorr.
R_cross = fftshift(R_cross);

figure;
theme(gcf, "light");
plot(lags, R_cross);
grid on;
xlabel('Lag');
ylabel('Cross-correlation');
title(sprintf("Cross-correlation of m-sequence with N = %d", L));
axis('tight');

%% A2.1.2.9 Compare the two results
% As we can see the cross correlation of both codes seems like a noisy waveform, in the first one (gold codes) we have the shape that is modeled in
% three levels, one at -t (-17) one at t - 2 (15) and one at -1.
% 
% In the other Cross correlation plot the shape of the wave is a lot more noisy than the previous one, in terms of levels we don't have the classic three
% levels of gold codes simply because the sequence are not gold codes, even
% without using the same function 

