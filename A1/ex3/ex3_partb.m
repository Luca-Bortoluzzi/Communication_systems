clc;
clear;
close all;

%% Exercise 3b

%Generate a non-systematic matrix generator G
G = [1 1 0 1 0;
    0 1 1 0 1;
    1 0 1 1 1];

rank = rank(G);

k = 3; % information bits length
n = 6; % codeword length (k + rank(G))

v = dec2bin(0:2^k-1) - '0';     % creating vectors minus zero ASCII

c = mod(v*G, 2);

y = dec2bin(0:2^n-1) - '0';     % possible received vectors

% we need to calculate the hamming distance to choose the right decision
% region

% hamming weight
Hw_c = sum(c, 2);
Hw_y = sum(y, 2);
% min hamming distance
%Hd = min(Hw(Hw>0));

Hd = zeros(64, 1);

for i=1:2^n
    Hd(i, 1) = min(diff(Hw_c(:), Hw_y(i, 1), 2))
end
