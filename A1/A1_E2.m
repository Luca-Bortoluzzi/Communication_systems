clear all
clc
close all

% Definisco il polinomio generatore g(D) = (D+1)*(D^3 + D + 1) = (D^4 + D^3 + D^2 + 1) e ne ricavo
% il vettore dei coefficienti
g_pol = [1 1 1 0 1];
r = length(g_pol) - 1; %grado del polinomio generatore

for k=2:4 
    m = k + r; %calcolo della lunghezza della codeword
    num_seq = 2^k;
    V = dec2bin(0:num_seq-1, k) - '0'; % Matrice 2^k x k

    C = zeros(num_seq, m);
    weights = zeros(num_seq, 1);
    
    fprintf('--------------------------------------------------\n');
    fprintf(' Risultati per k = %d (Lunghezza codice m = %d)\n', k, m);
    fprintf('--------------------------------------------------\n');
    
    if k == 2
        fprintf('v_vec\t\tc_vec\t\t\tWeight\n');
    end

end

%% struttura base CRC ESEMPIO (non c'entra nulla on la task)
clc
clear all
close all
% 1\. PARAMETRI DEL SISTEMA E POLINOMIO GENERATORE g(D) 
% Messaggio di informazione v di k = 4 bit 
v = [0 1 0 0]; 
k = length(v);

% Polinomi su F2: % p(D) = D^3 + D + 1 (Polinomio primitivo di grado 3: ) 
% (D + 1) = (Fattore per parità pari) 
% g(D) = (D + 1) * p(D) = D^4 + D^3 + D^2 + 1 -> Vettore: 
g = 'z^4 + z^3 + z^2 + 1';
crc_encoder = comm.CRCGenerator('Polynomial', g);
crc_detector = comm.CRCDetector('Polynomial', g); 

r = length(g) - 1;


% 2\. CODIFICA SISTEMATICA DEL CRC
% Passo A: Traslazione a sinistra di r posizioni -> v(D) * D^r 
v_shifted = [v, zeros(1, r)];
%Passo B: Divisione polinomiale su F2 tra v(D)*D^r e g(D) per trovare il
%resto r(D)
[~, resto] = gfdeconv(v_shifted, g);

%Passo C: Costrizione della parola di codive c = [v | resto]
c = [v, resto];
fprintf('--- TRASMETTITORE (CODIFICA) ---\n');
fprintf('Bit di parità calcolati r(D): [%s]\n', num2str(resto));
fprintf('Parola di codice trasmessa c: [%s] (Lunghezza n=%d)\n\n', num2str(c), k+r);

%% 
% =========================================================================
% Exercise 2 (Part 3) - CRC Burst Error Detection con comm.CRCGenerator
% =========================================================================

clear; clc;

% 1. Definizione del polinomio e creazione degli oggetti CRC
% Nota: usiamo gli apici dritti singoli (')
g = 'z^4 + z^3 + z^2 + 1';

crc_encoder = comm.CRCGenerator('Polynomial', g);
crc_detector = comm.CRCDetector('Polynomial', g);

% 2. Messaggio di informazione con k = 4 bit (vettore colonna)
v = [1; 0; 1; 1];

% 3. Codifica della parola di codice (lunghezza n = k + r = 4 + 4 = 8)
codeword = crc_encoder(v);

fprintf('--- TRASMETTITORE ---\n');
fprintf('Messaggio v:           [%s]^T\n', num2str(v'));
fprintf('Codeword c trasmessa:  [%s]^T (n = %d)\n\n', num2str(codeword'), length(codeword));

% =========================================================================
% 4. Definizione e Test dei Burst Errors (Punti 14-16)
% =========================================================================

% Errore e1 (Rilevato): Burst di lunghezza B = 3 <= deg(g) = 4
% Modifica i primi 3 bit
e1 = [1; 1; 1; 0; 0; 0; 0; 0];
rx1 = mod(codeword + e1, 2);
[~, err1] = crc_detector(rx1);

% Errore e2 (Non Rilevato): Multiplo esatto del polinomio generatore g(D)
% e2(D) = D^3 * g(D) = D^7 + D^6 + D^5 + D^3 -> [1 1 1 0 1 0 0 0]
e2 = [1; 1; 1; 0; 1; 0; 0; 0];
rx2 = mod(codeword + e2, 2);
[~, err2] = crc_detector(rx2);

% =========================================================================
% 5. Stampa dei Risultati
% =========================================================================

fprintf('--- RICEVITORE (DETECTOR) ---\n');
fprintf('Test 1 - Errore e1 (Burst B=3 <= 4):\n');
fprintf('  Vettore e1: [%s]^T\n', num2str(e1'));
fprintf('  Esito err:  %d (1 = Errore Rilevato) -> OK!\n\n', err1);

fprintf('Test 2 - Errore e2 (e2(D) = D^3 * g(D)):\n');
fprintf('  Vettore e2: [%s]^T\n', num2str(e2'));
fprintf('  Esito err:  %d (0 = Errore NON Rilevato / Sindrome Nulla) -> OK!\n', err2);


