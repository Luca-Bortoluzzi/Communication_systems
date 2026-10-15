clc;
clear;
% 3.1
I = eye(5);
A = [1 1 1 0 0;
    0 1 1 1 0;
    0 0 1 1 1;
    1 0 0 1 1;
    1 1 0 0 1;];
G1 = [I A];

P = A;
I = eye(5);

HT = [P; I];

% 3.2
m = [1 0 0 1 0];
y1 = m*G1;
% s1 = mod(y1*HT, 2);

e = [1 0 0 0 0 0 0 0 0 0];
y2 = mod(y1+e, 2);
% s2 = mod(y2*HT, 2);

% 3.3
HT_1 = HT;
HT_1(:, 1) = mod(HT(:, 1) + HT(:, 2), 2);
%s1 = mod(y1*HT_1, 2);
%s2 = mod(y2*HT_1, 2);

%3.4
%HT_2 = HT;
HT_2 = [HT, HT(:, 1), HT(:, 2)];

s1 = mod(y1*HT_2, 2);
s2 = mod(y2*HT_2, 2);
