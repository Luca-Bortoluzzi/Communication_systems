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

prbs_seq = zeros(1, N);     % Sequenza seriale di bit (1..1023)
states_history = zeros(N, m); % Storico di tutti gli stati

for k = 1:N
    % 1. Salva lo stato corrente e il bit in uscita (cella 10)
    states_history(k, :) = state;
    prbs_seq(k) = state(10);

    % 2. Calcola il bit di retroazione (XOR tra cella 10 e cella 7)
    feedback = xor(state(10), state(7));

    % 3. Scorrimento a destra: il feedback entra in Q1, Q10 esce
    state = [feedback, state(1:end-1)];
end

% Verifica dei risultati
fprintf('Primi 20 bit generati: ');
fprintf('%d', prbs_seq(1:20));
N0 = sum(prbs_seq == 0);
N1 = sum(prbs_seq == 1);
fprintf('\nNumero di 1 generati: %d (atteso: %d)\n', N1, 2^(m-1));
fprintf('Numero di 0 generati: %d (atteso: %d)\n', N0, 2^(m-1) - 1);

% Verifica che lo stato dopo 1023 cicli torni esattamente al seed iniziale
if isequal(state, states_history(1, :))
    disp('Verifica periodo massimo: SUPERATA (ritorno al seed dopo 1023 cicli).');
end

% Conversion to bipolar sequence (0 -> -1, 1 -> +1) and plot
bipolar_seq = 2 * prbs_seq - 1;     % zero and ones to -1 and +1
figure;
hold on;
stairs(1:20, bipolar_seq(1:20), 'LineWidth', 1); % square wave

grid on;
ylim([-1.5 1.5]);
xlim([0 21]);
xticks(1:20);
yticks([-1 0 1]);
xlabel('Bit index (n)');
ylabel('Bipolar Amplitude');
title('First 20 bit of the bipolar sequence LFSR (P(x) = x^{10} + x^7 + 1)');

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

% 1. Allinea la sequenza in modo che non inizi e finisca con lo stesso bit

seq_runs = prbs_seq();      % to not lose the original value when rotate with circshift
shift_amount = 0;

while seq_runs(1) == seq_runs(end)
    seq_runs = circshift(prbs_seq, 1);
    shift_amount = shift_amount + 1;
end

% 2. Conta la lunghezza di ogni run con un ciclo semplice
NR0 = zeros(1, m);   % zeros counter
NR1 = zeros(1, m);   % ones counter

run_lengths = [];
run_starts = [];
run_bits = [];

count = 1;
current_start = 1;

for k = 2:N
    if seq_runs(k) == seq_runs(k-1)
        count = count + 1;          % Il bit è uguale al precedente: allungo il run
    else
        run_lengths(end+1) = count;
        run_starts(end+1) = current_start;
        run_bits(end+1) = seq_runs(k-1);

        if seq_runs(k-1) == 0
            NR0(count) = NR0(count) + 1;
        else
            NR1(count) = NR1(count) + 1;
        end

        % Resetta per il nuovo run
        count = 1;
        current_start = k;
    end
end

% Salva l'ultimo run della sequenza
run_lengths(end+1) = count;
run_starts(end+1)  = current_start;
run_bits(end+1)    = seq_runs(N);

if seq_runs(N) == 0
    NR0(count) = NR0(count) + 1;
else
    NR1(count) = NR1(count) + 1;
end

% 3. Costruzione e stampa della tabella
NR_i = NR0 + NR1;
NR_T = sum(NR_i);
ideal_vals = (NR_T ./ (2 .^ (1:m)))';

fprintf('Numero totale di run NR(T) = %d\n\n', NR_T);
T_runs = table((1:m)', NR0', NR1', NR_i', ideal_vals, ...
    'VariableNames', {'Length_i', 'NR0_i', 'NR1_i', 'NR_i', 'Ideal_NR_over_2i'});
disp(T_runs);

% 4. Grafico: Run length vs. Starting point (colori diversi per 0 e 1)
idx_0 = (run_bits == 0);
idx_1 = (run_bits == 1);

figure;
stem(run_starts(idx_0), run_lengths(idx_0), 'cyan', 'filled', 'DisplayName', 'Run di 0');
hold on;
stem(run_starts(idx_1), run_lengths(idx_1), 'o', 'filled', 'DisplayName', 'Run di 1');

grid on;
xlim([1, N]);
ylim([0, m + 1]);
yticks(1:m);
xlabel('Indice di inizio nella sequenza (Starting point)');
ylabel('Lunghezza del run (i)');
title('Lunghezza dei run vs. Punto di inizio nella sequenza (m = 10)');
legend('Location', 'best');


%% --- Periodic autocorrelation ---
% mean value E[ai] = 0
mean_val = mean(bipolar_seq);
fprintf('Media simboli bipolari E[a_i]: %.5f (\x2248 0)\n', mean_val);

% Intervallo di tau richiesto dalla consegna: -(N-1)/2 <= tau <= +(N-1)/2
max_tau = (N - 1) / 2;            % Per N = 1023 vale 511
lags    = -max_tau : max_tau;     % Vettore dei ritardi da -511 a +511
R_tau   = zeros(size(lags));

% Ciclo for che fa scorrere tau
for k = 1:length(lags)
    tau = lags(k);
    shifted_seq = circshift(bipolar_seq, tau);      % Calcola v'(i - tau) periodico
    R_tau(k)    = sum(bipolar_seq .* shifted_seq);  % Somma da i=0 a N-1
end

% Plot dell'autocorrelazione
figure;
plot(lags, R_tau, 'LineWidth', 1.5);
grid on;
xlim([-max_tau, max_tau]);
ylim([-100, N + 100]);
xlabel('\tau');
ylabel('R_{v''}(\tau)');
title('Autocorrelazione Periodica della sequenza bipolare (m = 10, N = 1023)');