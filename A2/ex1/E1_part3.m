clc;
clear;
close all;

N = 127;
K_vals = 1:200;
num_K = length(K_vals);

% Pre-allocazione dei vettori per i quattro bound
welch_gen    = zeros(1, num_K); % Welch generale (s <= 10)
welch_s1     = zeros(1, num_K); % Welch caso particolare (s = 1)
welch_approx = zeros(1, num_K); % Welch approssimazione (sqrt(N))
sidelnikov   = zeros(1, num_K); % Sidelnikov bound

% Intervalli ammissibili per s
s_welch = 1:10;
max_s_sid = min(10, ceil(2*N/5) - 1);
s_sidelnikov = 0:max_s_sid;

for i = 1:num_K
    K = K_vals(i);

    % 3. Approssimazione Welch (indipendente da K)
    welch_approx(i) = ceil(sqrt(N));

    % --- 1. Welch generic (s <= 10) ---
    temp_w = zeros(size(s_welch));
    for idx = 1:length(s_welch)
        s = s_welch(idx);
        binom = nchoosek(N + s - 1, s);
        arg = (1 / (K*N - 1)) * ((K*N) / binom - 1);
        arg = max(0, arg);
        temp_w(idx) = N * (arg ^ (1 / (2*s)));
    end
    welch_gen(i) = ceil(max(temp_w));

    % --- 2. Welch (s = 1) ---
    welch_s1(i) = ceil(N * sqrt((K - 1) / (K*N - 1)));

    % --- 4. Sidelnikov bound (0 <= s <= 10, s < 2N/5) ---
    temp_sid = zeros(size(s_sidelnikov));
    for idx = 1:length(s_sidelnikov)
        s = s_sidelnikov(idx);
        term1 = (2*s + 1) * (N - s);
        term2 = (s * (s + 1)) / 2;
        term3 = (2^s * N^(2*s + 1)) / (K * factorial(2*s) * nchoosek(N, s));

        inner_arg = max(0, ceil(term1 + term2 - term3));
        temp_sid(idx) = ceil(sqrt(inner_arg));
    end
    sidelnikov(i) = max(temp_sid);
end

m = 7;              % numero di celle (LFSR)
N = 2^m - 1;        % lunghezza della singola sequenza (127)
K_gold = N + 2;          % numero totale di sequenze Gold (129)

% 1. Generazione dell'intero set di K = 129 sequenze bipolari (ciascuna lunga N = 127)
X_bipolar = zeros(K_gold, N);
idx_values = -2:(N - 1); % Da -2 a 126 -> 129 sequenze

for k = 1:K_gold
    goldseq = comm.GoldSequence('FirstPolynomial', 'x^7+x^3+1', ...
        'SecondPolynomial', 'x^7+x+1', ...
        'FirstInitialConditions',  [0 0 0 0 0 0 1], ...
        'SecondInitialConditions', [0 0 0 0 0 0 1], ...
        'Index', idx_values(k), ...
        'SamplesPerFrame', N);

    x = goldseq()';
    X_bipolar(k, :) = 1 - 2*x; % Conversione in {+1, -1}
end

% 2. Calcolo di r_M tramite FFT su tutte le coppie (auto e cross-correlazione)
X_fft = fft(X_bipolar, [], 2);
rM_gold = 0;

for i = 1:K_gold
    for j = i:K_gold
        R_periodic = round(real(ifft(X_fft(i, :) .* conj(X_fft(j, :)))));

        if i == j
            % Autocorrelazione: escludiamo il picco a lag 0 (che vale N = 127)
            R_periodic(1) = 0; 
        end

        rM_gold = max(rM_gold, max(abs(R_periodic)));
    end
end

disp(['r_M calcolato = ', num2str(rM_gold)]); % Restituirà 17


% Grafico comparativo
figure;
theme(gcf, "light");
plot(K_vals, sidelnikov,   'LineWidth', 2.2); hold on;
plot(K_vals, welch_gen,    'LineWidth', 2);
plot(K_vals, welch_s1,     '--', 'LineWidth', 2);
plot(K_vals, welch_approx, ':',  'LineWidth', 2);

% Aggiunta del singolo punto (129, 17)
plot(K_gold, rM_gold, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

grid on;
xlabel('Number of sequences (K)');
ylabel('Lower bound on r_M');
title(['Welch vs Sidelnikov Bounds compare (N = ', num2str(N), ')']);
legend('Sidelnikov bound', ...
    'Welch generic(s \leq 10)', ...
    'Welch (s = 1)', ...
    'Welch approx (\lceil\surdN\rceil)', ...
    ['Gold sequences (K = ', num2str(K_gold), ', r_M = ', num2str(rM_gold), ')'], ...
    'Location', 'southeast');

% Explanation of how r_M was computed:
% We generated the entire set of K = 129 Gold sequences of length N = 127
% by varying the 'Index' parameter from -2 to 126 (N - 1). Specifically,
% Index = -2 and Index = -1 correspond to the two base m-sequences generated
% by the preferred pair of primitive polynomials x^7 + x^3 + 1 and x^7 + x + 1,
% while Index = 0 to 126 correspond to their 127 cyclic-shifted combinations.
% After converting the sequences to bipolar form (+1, -1), we computed the
% periodic autocorrelation (for non-zero lags, excluding the peak of N = 127
% at lag 0) and periodic cross-correlation across all sequence pairs and lags.
% Taking the maximum absolute value yields r_M = 17, which was then plotted
% as the point (129, 17) on the bound figure.

%% --- comment the results ---

% 1. Welch bounds: For K <= N (K <= 127), the general Welch bound (s <= 10)
%    and its simplified version (s = 1) coincide and asymptotically approach
%    the constant approximation ceil(sqrt(N)) = 12. However, for K > N, the
%    s = 1 bound saturates at 12, whereas the general Welch bound keeps growing
%    because higher-order terms (s >= 2) become dominant.
%
% 2. Sidelnikov bound: The Sidelnikov bound lies strictly above all three
%    Welch curves across the entire range of K, confirming that it provides a
%    tighter lower bound by exploiting the binary bipolar nature of the sequences.
%    Due to the integer ceiling operations and the maximization over s, it remains
%    flat over intervals of K and increases in discrete steps.
%
% 3. Gold sequences optimality: For K = 129, the Gold sequence set achieves
%    r_M = 17, which sits just 1 unit above the Sidelnikov lower bound (16) and
%    well above the Welch bounds (13 and 12). This demonstrates that Gold
%    sequences are near-optimal with respect to the Sidelnikov bound while
%    offering a large family size (K = N + 2) with bounded three-valued correlation.