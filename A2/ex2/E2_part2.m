%% Exercise 2 part 2

clc;
clear;
close all;

Eb = 1;     % energy per bit (pulse unitary energy)

errors = 100;
EbN0_db = 0:8;     % Eb / N0
EbN0_lin = 10.^(EbN0_db / 10);  % Linear Conversion
P_a = 1/2 * erfc(sqrt(EbN0_lin));

v1 = rand([0, 1])
c1 = dec2bin()