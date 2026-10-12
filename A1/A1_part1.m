clc;
clear;
close all;
I = eye(5);
A = [1 1 1 0 0;
    0 1 1 1 0;
    0 0 1 1 1;
    1 0 0 1 1;
    1 1 0 0 1;];
B = [1 1 0 0 0;
    0 1 1 1 0;
    0 0 1 1 1;
    1 0 0 1 1;
    1 1 0 0 1;];
G1 = [I A];
G2 = [I B];
k = 5;
%% 1 
v = dec2bin(0:2^k-1) - '0';
%% 2
c1 = mod(v*G1, 2);
c2 = mod(v*G2, 2);
%% 3 and 4
%variable define
hw1=zeros(32, 1);
hw2=zeros(32, 1);
mult1=zeros(11,1);
mult2=zeros(11,1);
for i=1:32
    %hamming weight
    hw1(i, 1) = sum(c1(i, : ));
    hw2(i, 1) = sum(c2(i, :));
    %Multiplicities
    a1 = hw1(i, 1); %%temp variables
    a2 = hw2(i, 1);
    mult1(a1+1, 1) = mult1(a1+1, 1) + 1;
    mult2(a2+1, 1) = mult2(a2+1, 1) + 1;
end
%% 5
% we don't calculate the hamming distance of vector row 32 because it's
% redundant since it is calculated by all the precedent vectors so on the
% 32x32 matrix hd it would all be zeroes
min1 = zeros(31,1);
min2 = zeros(31,1);
hd1 = zeros(32, 32);
hd2 = zeros(32, 32);
for i=1:31
    for j= i+1:32
        hd1(i,j) = sum(xor(c1(i,1:10), c1(j,1:10))) ;
        hd2(i,j) = sum(xor(c2(i,1:10), c2(j,1:10))) ;

    end
    min1(i,:) = min(setdiff(unique(hd1(i,:)), 0));
    min2(i,:) = min(setdiff(unique(hd2(i,:)), 0));
end



