clc;
clear;

I = eye(5);

P1 = [1 1 1 0 0;
    0 1 1 1 0;
    0 0 1 1 1;
    1 0 0 1 1;
    1 1 0 0 1];

G1 = [I, P1];

P2 = [1 1 0 0 0;
    0 1 1 1 0;
    0 0 1 1 1;
    1 0 0 1 1;
    1 1 0 0 1];

G2 = [I, P2];

v = dec2bin(0:(2^5 - 1)) - '0';

c1 = mod(v*G1, 2);
c2 = mod(v*G2, 2);
hw1 = zeros(32, 1);
hw2 = zeros(32, 1);
mult1 = zeros(11, 1);
mult2 = zeros(11, 1);

for i=1:32
    hw1(i,1) = nnz(c1(i, :));
    hw2(i,1) = nnz(c2(i, :));
end
for j=1:11
    mult1(j,1) = mult1(j, 1) + (sum(hw1(:,1) == (j-1)));
    mult2(j,1) = mult2(j, 1) + (sum(hw2(:,1) == (j-1)));
end

d1 = zeros(31, 32);
d2 = zeros(31, 32);

for i=1:31
    for k=i+1:32
        d1(i,k) = sum(xor(c1(i, :), c1(k, :)));
        d2(i,k) = sum(xor(c2(i, :), c2(k, :)));
    end
end

d_min1 = zeros(31,1);
d_min2 = zeros(31,1);
for i=1:31
    d_min1(i) = min(d1(i, i+1:end));
    d_min2(i) = min(d2(i, i+1:end));
end

min1 = min(d_min1);
min2 = min(d_min2);

mult_min1 = 0;
mult_min2 = 0;

for j=1:31
    mult_min1 = mult_min1 + (sum(d_min1(j) == (min1)));
    mult_min2 = mult_min2 + (sum(d_min2(j) == (min2)));
end

