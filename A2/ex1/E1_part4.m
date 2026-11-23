%1.4
clc;
clear;
close all;
m=7;
K_gold=2^m+1;
N_originale=2^m-1;

% Generazione dell'intero set di K = 129 sequenze bipolari (ciascuna lunga N = 127)
X_bipolar_127 = zeros(K_gold, N_originale);
idx_values = -2:(N_originale - 1); % Da -2 a 126 -> 129 sequenze

for k = 1:K_gold
    goldseq = comm.GoldSequence('FirstPolynomial', 'x^7+x^4+1', ...
        'SecondPolynomial', 'x^7+x+1', ...
        'FirstInitialConditions',  [0 0 0 0 0 0 1], ...
        'SecondInitialConditions', [0 0 0 0 0 0 1], ...
        'Index', idx_values(k), ...
        'SamplesPerFrame', N_originale);

    x = goldseq()';
    X_bipolar_127(k, :) = 1 - 2*x; % Conversione in {+1, -1}
end

%1.4.1 Rimozione ultimo chip e calcolo di R_m
X_bipolar = X_bipolar_127(:,1:end-1);
N=size(X_bipolar,2);
R_A=0;
R_C=0;
for i= 1:K_gold
    u= X_bipolar(i,:);
    auto_corr = round(ifft(fft(u) .* conj(fft(u))));
    R_A = max(R_A, max(abs(auto_corr(2:end))));
    for j=i+1:K_gold
        v=X_bipolar(j,:);
        cross_corr = round(ifft(fft(u) .* conj(fft(v))));
        R_C = max(R_C, max(abs(cross_corr)));
    end
end
R_M=max(R_A,R_C);

K_vec=1:200;
welch_sqrt = ceil(sqrt(N)) * ones(size(K_vec));
welch_s1   = zeros(size(K_vec));
welch_orig = zeros(size(K_vec));
sidelnikov = zeros(size(K_vec));
for idx=1:length(K_vec)
    K_curr=K_vec(idx);
    % Welch formula semplificata (s = 1)
    arg_s1 = (K_curr - 1) / (K_curr * N - 1);
    welch_s1(idx) = ceil(N * sqrt(max(0, arg_s1)));

    % Welch bound originale (max su 1 <= s <= 10)
    max_w = 0;
    for s = 1:10
        denom = nchoosek(N + s - 1, s);
        arg = (1 / (K_curr * N - 1)) * (K_curr * N / denom - 1);
        if arg > 0
            max_w = max(max_w, N * (arg)^(1 / (2 * s)));
        end
    end
    welch_orig(idx)=ceil(max_w);
    max_sid = 0;
    s_max = min(10, ceil(2 * N / 5) - 1);
    for s = 0:s_max
        log_num = s * log(2) + (2*s + 1) * log(N);
        log_den = log(K_curr) + gammaln(2*s + 1) + gammaln(N + 1) - gammaln(s + 1) - gammaln(N - s + 1);
        term = exp(log_num - log_den);
        arg = (2*s + 1)*(N - s) + s*(s + 1)/2 - term;
        if arg > 0
            max_sid = max(max_sid, sqrt(arg));
        end
    end
    sidelnikov(idx) = ceil(max_sid);
end

%1.4.2
figure;
plot(K_vec, sidelnikov, 'r-', 'LineWidth', 1.2); hold on;
plot(K_vec, welch_orig, 'b-o', 'LineWidth', 1);
plot(K_vec, welch_s1, 'c-', 'LineWidth', 1.2);
plot(K_vec, welch_sqrt, 'k-', 'LineWidth', 1);
grid on,axis('tight');
xlabel('K');
ylabel('bound'),ylim([0 22]);
title('Welch and Sidelnikov bounds for N = 126');
legend('Sidelnikov binary', 'Welch original', 'Welch formula', 'Welch sqrt(N)', 'Location', 'northwest');
exportgraphics(gcf, 'A2_1_15.png', 'Resolution', 300);

%1.4.3
figure;
plot(K_vec, sidelnikov, 'r-', 'LineWidth', 1.2); hold on;
plot(K_vec, welch_orig, 'b-o', 'LineWidth', 1);
plot(K_vec, welch_s1, 'c-', 'LineWidth', 1.2);
plot(K_vec, welch_sqrt, 'k-', 'LineWidth', 1);
plot(K_gold, R_M, 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
grid on,axis('tight');
xlabel('K');
ylabel('bound');
title('Welch and Sidelnikov bounds for N = 126');
legend('Sidelnikov binary', 'Welch original', 'Welch formula', 'Welch sqrt(N)' ...
    ,sprintf('R_M = %d (K = %d)', R_M, K_gold), 'Location', 'northwest');
ylim([0 max(R_M + 4, 22)]);
exportgraphics(gcf, 'A2_1_16.png', 'Resolution', 300);

%1.4.4

% Dal confronto con i bound teorici calcolati per N = 126 e K = 129, si
% osserva che il valore massimo di correlazione ciclica r_M ottenuto
% (r_M = 26) risulta sensibilmente superiore rispetto al bound binario di
% Sidelnikov (pari a 16) e ai bound di Welch (pari a 12).
%
% A differenza del caso originale con N = 127 — in cui le sequenze di Gold
% raggiungono l'ottimalita' rispetto al bound di Sidelnikov (r_M = 17)
% grazie alla struttura algebrica dei campi di Galois GF(2^m) e alla
% scelta di una preferred pair di polinomi primitivi — la rimozione anche
% di un solo chip distrugge la periodicita' nativa N = 2^m - 1 e il
% bilanciamento tra +1 e -1 delle sequenze.
%
% Di conseguenza, si perde completamente la proprieta' di correlazione
% ciclica a tre valori {-t, -1, t-2} e i picchi di autocorrelazione fuori
% fase e cross-correlazione aumentano drasticamente, rendendo le sequenze
% troncate non ottimali per applicazioni CDMA asincrone.

