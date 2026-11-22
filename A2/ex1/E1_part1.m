%% Exercise 1 Assignment 2

close all;
clear;
clc;
%% A2.1.1

% code from the slides
m = 7;      % number of cells (FLSR)
L = 2^m - 1;
n = 0:L-1;

goldseq = comm.GoldSequence('FirstPolynomial','x^7+x^4+1',...
    'SecondPolynomial','x^7+x+1',...
    'FirstInitialConditions',[0 0 0 0 0 0 1],...
    'SecondInitialConditions',[0 0 0 0 0 0 1],...
    'Index',1,'SamplesPerFrame',L);

x0 = goldseq()';
x0b = 1 - 2*x0;


R0 = ifft(fft(x0b).*conj(fft(x0b))); %non-normalized periodic autocorr.
R0 = fftshift(R0);

%% A2.1.2
figure;
theme(gcf, "light");
title("Correlation of Gold code with N = 127")
plot(n,R0);
axis('tight');


%% A2.1.3 repeat for i = -1 and i = -2

release(goldseq);        % Sblocca l'oggetto
goldseq.Index = -1;      % Cambia solo l'indice
x1 = goldseq()';       % Genera la nuova sequenza
x1b = 1 - 2*x1;

R1 = ifft(fft(x1b).*conj(fft(x1b))); %non-normalized periodic autocorr.
R1 = fftshift(R1);

figure;
theme(gcf, "light");
plot(n,R1);
title("Correlation of Gold code with N = 127")
axis('tight');


release(goldseq);        % Sblocca l'oggetto
goldseq.Index = -1;      % Cambia solo l'indice
x2 = goldseq()';       % Genera la nuova sequenza
x2b = 1 - 2*x2;


R2 = ifft(fft(x2b).*conj(fft(x2b))); %non-normalized periodic autocorr.
R2 = fftshift(R2);

figure;
theme(gcf, "light");
title("Correlation of Gold code with N = 127")
plot(n,R2);
axis('tight');

%% comment the results A2.1.4
% come si può vedere dai plot, la funzione di autocorrelazione utilizzando come indice i = 1, è un grafico su 3 livelli eccetto per l'origine che ha correlazione
% massima pari a L = 2^m - 1 = 127. I tre livelli sono a -17, -1 e 15 sull'asse delle ordinate.
% Un'osservazione importante è la minima quasi nulla differenza tra le funzioni di autocorrelazione computate ad indici -1 e -2; 
% i grafici appaiono identici, ma ci sono delle differenze tra i codici bit per bit, valore pari a -1 per qualsiasi valore di N escluso il valore
% lag = 0 ove il grafico vale 127 come visto precedentemente.




