% =========================================================================
% Exercise 2 – CRC Codes and Weight Distribution
% =========================================================================
clc;        % Pulisce la finestra di comando
clear;      % Rimuove tutte le variabili dalla memoria per evitare conflitti
close all;  % Chiude eventuali figure aperte



% --- 1. & 2. Definizione del Polinomio Generatore g(z) = (z + 1) * p(z) ---
% Il polinomio base è p(z) = z^3 + z + 1, che è un polinomio primitivo di grado r_p = 3.
% Moltiplicando per (z + 1) otteniamo il polinomio generatore completo:
% g(z) = (z + 1)*(z^3 + z + 1) = z^4 + z^3 + z^2 + 1 (di grado r = 4).
g_str = 'z^4 + z^3 + z^2 + 1';

% Inizializzazione del CRC encoder del Communications Toolbox per la
% generazione del CRC
crc_encoder = comm.CRCGenerator('Polynomial', g_str);



% --- 3. Parametri iniziali del codice e sequenze d'ingresso ---
k = 2;                             % Numero di bit di informazione da trasmettere
r = 4;                             % Grado del polinomio generatore g(z) (numero di bit di parità)
n = k + r;                         % Lunghezza totale della parola di codice (n = 2 + 4 = 6)

% Generazione di tutte le 2^k possibili sequenze d'ingresso v come vettori riga binari.
v_matrix = dec2bin(0:2^k-1) - '0';



% --- 4. Costruzione della matrice generatrice sistematica G e codifica ---
% Inizializzazione della matrice G di dimensione (k x n)
G = zeros(k, n);
I_k = eye(k);                      % Matrice identità di dimensione k x k

% Le righe di G sono ricavate passando le colonne dell'identità I_k nel sistema CRC
for i = 1:k
    v_col = I_k(i, :)';            % Estrazione dell'i-esima riga dell'identità e trasposizione
    codeword_col = crc_encoder(v_col); % Il CRC appende i bit di parità in coda al vettore
    G(i, :) = codeword_col';       % Salvataggio della parola di codice ottenuta come riga di G
end

disp('=== MATRICE GENERATRICE SISTEMATICA G ===');
disp(G);

% Calcolo di tutte le 2^k codewords tramite il prodotto matriciale v * G in modulo 2
c_matrix = mod(v_matrix * G, 2);

% Calcolo del peso di Hamming w_H per ogni parola di codice
w_H = sum(c_matrix, 2);



% --- 5. Stampa dei dati per la tabella preliminare ---
disp(' ');
disp('=== PUNTO 5: VALORI PER LA TABELLA ===');
disp('Vettori v (informazione):');
disp(v_matrix);
disp('Codewords c (parole di codice):');
disp(c_matrix);
disp('Pesi di Hamming w_H per ciascuna parola:');
disp(w_H);



% --- 6. Calcolo della distribuzione dei pesi A_i ---
% Inizializzazione del vettore A per contare quante codewords hanno peso i (con i da 0 a n)
A = zeros(1, n + 1);
for i = 0:n
    % L'indice nell'array parte da 1, quindi A(i+1) memorizza il conteggio per peso = i
    A(i+1) = sum(w_H == i);
end

disp('=== PUNTO 6: DISTRIBUZIONE DEI PESI A_i ===');
for i = 0:n
    fprintf('A_%d = %d\n', i, A(i+1));
end



% --- 7. Commento sui risultati ottenuti per k = 2 ---
% La distanza minima d_min è il peso minimo tra tutte le codewords non nulle
d_min = min(w_H(w_H > 0));

fprintf('\n=== PUNTO 7: COMMENTO SUI RISULTATI ===\n');
fprintf('- Peso minimo del codice (d_min): %d\n', d_min);
fprintf('- Tutti i pesi di Hamming sono numeri PARI (w_H in {0, 4, 6}).\n');
fprintf('- La presenza del fattore (z+1) nel polinomio generatore assicura che\n');
fprintf('  ogni parola di codice abbia un numero pari di 1, consentendo di\n');
fprintf('  rilevare qualsiasi errore di canale con peso DISPARI.\n');



% --- Punto 8: Ricerca iterativa aumentando k finché non compare w_H = 2 ---
% flag booleano per controllare il ciclo while
trovato_peso_2 = false;
k_curr = k;                        % Inizio della ricerca a partire da k = 2

% Strutture dati adibite al salvataggio
risultati_k = [];                  % Salva i valori di k analizzati
risultati_A = cell(0);             % Cell array per salvare i vettori A_i di lunghezza variabile

while ~trovato_peso_2
    n_curr = k_curr + r;           % Nuova lunghezza del codice: n = k + 4
    
    % Generazione di tutte le 2^k combinazioni di informazione per il k corrente
    v_curr = dec2bin(0 : 2^k_curr - 1) - '0';
    
    % Costruzuione della matrice G aggiornata per la nuova dimensione k_curr x n_curr
    G_curr = zeros(k_curr, n_curr);
    I_curr = eye(k_curr);
    for i = 1:k_curr
        v_col = I_curr(i, :)';
        cw_col = crc_encoder(v_col);
        G_curr(i, :) = cw_col';
    end
    
    % Calcolo delle codewords e i relativi pesi di Hamming in GF(2)
    c_curr = mod(v_curr * G_curr, 2);
    w_H_curr = sum(c_curr, 2);     % Somma lungo le righe per contare i bit a 1
    
    % Calcolo dell vettore di distribuzione dei pesi A_i per questa configurazione
    A_curr = zeros(1, n_curr + 1);
    for i = 0:n_curr
        A_curr(i+1) = sum(w_H_curr == i);
    end
    
    % Salvataggio nei risultati dei valori trovati in questo step
    risultati_k = [risultati_k, k_curr];
    risultati_A{end+1} = A_curr;
    
    % Verifica dell'esistenza di codewords con peso pari a 2.
    % Poiché A_curr(1) corrisponde a peso 0 e A_curr(2) a peso 1,
    % il valore di A_2 (cioè il numero di codewords di peso 2) si trova in A_curr(3).
    if A_curr(3) > 0
        trovato_peso_2 = true;     % Condizione di arresto: trovato w_H = 2!
    else
        k_curr = k_curr + 1;       % Altrimenti incremento k e continuo a cercare
    end
end



% --- Punto 9: Guida e valori calcolati per compilare la tabella a mano ---
fprintf('\n=========================================================\n');
fprintf('=== PUNTO 9: GUIDA PER COMPILARE LA TABELLA A MANO ===\n');
fprintf('=========================================================\n');
fprintf('Per tracciare la tabella della distribuzione dei pesi su carta, organizzare\n');
fprintf('le colonne nel formato: [ k | n | A_0 | A_1 | A_2 | ... | A_n ]\n\n');
fprintf('Riportare i seguenti dati calcolati dal codice:\n\n');



% ---------------------------------------------------------------------
% Ciclo di stampa pulita per trascrivere facilmente i dati sul foglio
% DA CANCELLARE POI
for idx = 1:length(risultati_k)
    k_val = risultati_k(idx);
    n_val = k_val + r;
    A_val = risultati_A{idx};
    
    fprintf('• Per k = %d (con n = %d):\n  ', k_val, n_val);
    for i = 0:n_val
        fprintf('A_%d = %d', i, A_val(i+1));
        if i < n_val
            fprintf(', ');
        end
    end
    fprintf('\n\n');
end
% ---------------------------------------------------------------------



% --- Punto 10: Dimostrazione e spiegazione teorico-matematica ---
fprintf('=========================================================\n');
fprintf('=== PUNTO 10: SPIEGAZIONE MATEMATICA (per la relazione) ===\n');
fprintf('=========================================================\n');
fprintf('Le codeword di peso 2 compaiono per la prima volta quando k = %d (ovvero n = %d).\n\n', k_curr, k_curr + r);
fprintf('Dimostrazione teorica step-by-step:\n');
fprintf('1. Una parola di codice c(z) ha peso di Hamming pari a 2 se contiene solo due bit a 1.\n');
fprintf('   In forma polinomiale si scrive come: c(z) = z^i + z^j = z^j * (z^(i-j) + 1), con 0 <= j < i < n.\n\n');
fprintf('2. Per le proprietà dei codici ciclici/CRC, c(z) è una codeword valida se e solo se\n');
fprintf('   è divisibile per il polinomio generatore g(z):\n');
fprintf('   g(z) | c(z)  ==>  (z + 1)*p(z) | z^j * (z^(i-j) + 1).\n\n');
fprintf('3. Poiché il fattore z^j è coprimo con g(z), la condizione di divisibilità si riduce a:\n');
fprintf('   g(z) deve dividere il polinomio (z^L + 1), dove L = i - j è la distanza tra i due bit a 1.\n\n');
fprintf('4. Il polinomio p(z) = z^3 + z + 1 è un polinomio PRIMITIVO di grado r_p = 3.\n');
fprintf('   Il suo periodo (o ordine) è pari a N = 2^(r_p) - 1 = 2^3 - 1 = 7.\n');
fprintf('   Per teoria dei campi finiti, p(z) divide (z^L + 1) se e solo se L è un multiplo di N = 7.\n\n');
fprintf('5. La minima distanza non nulla tra i due bit affinché p(z) divida (z^L + 1) è dunque L = 7.\n');
fprintf('   Il fattore (z + 1) divide sempre (z^7 + 1) in aritmetica modulo 2, quindi g(z) divide (z^7 + 1).\n\n');
fprintf('6. Per ospitare due bit a 1 distanti L = 7 (ovvero per avere il polinomio z^7 + 1),\n');
fprintf('   il grado della codeword deve essere almeno 7, il che richiede una lunghezza\n');
fprintf('   minima del codice pari a n = L + 1 = 7 + 1 = 8.\n\n');
fprintf('7. Sapendo che la lunghezza del codice è n = k + r con r = 4 (grado di g(z)):\n');
fprintf('   n >= 8  ==>  k + 4 >= 8  ==>  k >= 4.\n');
fprintf('   Ecco spiegato matematicamente perché la prima parola di codice di peso 2\n');
fprintf('   è comparsa esattamente per k = 4!\n');



% --- Punto 12: Considernado un CRC code con generatore ---
% polinomiale g(D) di grado r con termine noto g(0) = 1.
% Determinare la lunghezza massima di burst L che è sicuramente riconosciuta da questo CRC code.
%
% Un CRC con 'r' bit di ridondanza rileva con certezza assoluta qualsiasi
% errore a raffica di lunghezza L <= r.



% --- Punto 13: Dimostrazione matematica ---

% Un errore di lunghezza L viene espresso matematicamente come:
% e(D) = D^j * b(D)
% - b(D) è il Polinomio Base della raffica che definisce la composizione dei
% bit errati, con grado pari a degb(D)=L-1
% - D^j rappresenta lo scorrimento di 'j' posizioni all'interno del blocco
% trasmesso.
% un errore e(D) non viene
% rilevato dal ricevitore se e solo se è un multiplo esatto del polinomio
% generatore g(D).
% Poiché g(0) = 1, il termine D^j è coprimo (ossia non condivide fattore) con g(D); quindi l'errore
% è non rilevato se e solo se g(D) divide b(D).
% Se L <= r, allora deg(b(D)) = L - 1 <= r - 1 < r = deg(g), perciò b(D) non può
% essere un multiplo non nullo di g(D): ogni burst di lunghezza L <= r è rilevato.




% --- Punto 14: due pattern di errori a raffica tali che:
% - e1: è rilevato dal CRC con sindrome diversa da 0
% - e2: NON è rilevato dal CRC (sindrome nulla); ---

% vettori errore in colonna per compatibilità con CRCGenerator
% lunghezza burst L = 4 <= r  ---> rilevato
e1 = [1,0,1,1,0,0,0,0]';        % e1(D) = D^3 + D + 1
% lunghezza burst L = 5 = r+1 > r ---> NON rilevato
e2 = [1,1,1,0,1,0,0,0]';        % e2(D) = D^4 + D^3 + D^2 + 1 = g(D)




% --- Punto 15: motivazione dei due vettori di errore ---
% Per e1: 
%  il grado di deg(e1)=3  è minore del grado del polinomio generatore
%  deg(g)=4, di conseguenza la divisione e1(D)/g(D) dà un resto non nullo:
%       rem(e1(D,g(D)) = D^3 + D + 1 != 0
%  La sindrome è diversa da zero, quindi il CRC Rileva l'errore
% Per e2:
% e2(D) coincide con il polinomio generatore g(D), quindi è divisibile
% per g(D) e produce una sindrome nulla: il CRC non rileva l'errore.



% --- Punto 16: Scelta di una codeword e dimostrazione dei risultati ---
% Detector CRC inizializzato con lo stesso polinomio generatore g(D)
crc_det = comm.CRCDetector('Polynomial', g_str);

% sequenza di informazione a 4 bit (k = 4, dato che e1 ed e2 sono di lunghezza n = 8)
v_test = [1 0 1 1]';

% codeword c generata tramite l'encoder CRC
c = crc_encoder(v_test);

% Generazione dei segnali ricevuti sul canale (r = c + e mod 2)
r1 = mod(c + e1, 2);
r2 = mod(c + e2, 2);

% Verifica con il detector CRC
[~, err1] = crc_det(r1); % err1 = 1 se rilevato, 0 se non rilevato
[~, err2] = crc_det(r2);

% Stampa dei risultati
fprintf('=========================================================\n');
fprintf('=== PUNTO 16: RISULTATI VERIFICA BURST ERRORS ===\n');
fprintf('=========================================================\n');

fprintf('• Analisi Errore e1 (Burst L = 4):\n');
if err1 == 1
    fprintf('  -> ESITO: Errore RILEVATO dal CRC (Sindrome != 0).\n');
else
    fprintf('  -> ESITO: Errore NON rilevato.\n');
end

fprintf('\n• Analisi Errore e2 (Burst L = 5):\n');
if err2 == 1
    fprintf('  -> ESITO: Errore RILEVATO dal CRC.\n');
else
    fprintf('  -> ESITO: Errore NON RILEVATO dal CRC (Sindrome == 0).\n');
end