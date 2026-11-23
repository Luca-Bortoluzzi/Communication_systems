clc;
clear;
close all;

%% -- Exercise 2 - CDMA ---
%% --- Part 1: 2-PAM ---
% Plot the analytical BER vs Eb/N0 from 0 to 12 db

range = 0:12;
Eb = 1;     % energy per bit (pulse unitary energy)

N0 = randn();       % AWGN noise
args = Eb / N0;
P_a = 1/2 * erfc(sqrt(args));

figure;
theme(gcf, "light");
loglog(range, P_a);
