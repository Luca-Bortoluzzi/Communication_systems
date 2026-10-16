
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
%% 1.1
v = dec2bin(0:2^k-1) - '0';
%% 1.2
n = 10;                    % codeword lenght (10)
c1 = mod(v*G1, 2);
c2 = mod(v*G2, 2);
%% 1.3 and 1.4
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
    a2 = hw2(i, 1);         %ATTENTION hw2!!! not hw1!!!
    mult1(a1+1, 1) = mult1(a1+1, 1) + 1;
    mult2(a2+1, 1) = mult2(a2+1, 1) + 1;
end
%% 1.5
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

%% Part 2

%% 2.1 
% Evaluate PUE(p) for both codes using their weight spectra.
p_bsc = 0.5;     % Bernoulli probability for BSC

PUE1_bsc = 0;
PUE2_bsc = 0;

for w = 1:n
    PUE1_bsc = PUE1_bsc + mult1(w+1) * p_bsc^w * (1-p_bsc)^(n-w);
    PUE2_bsc = PUE2_bsc + mult2(w+1) * p_bsc^w * (1-p_bsc)^(n-w);
end

%% 2.2 
% Plot both curves on the same log–log axes for p ∈[10−6,10−1].
num = 1000;   %number of values for probability interval
p_interval = logspace(-6, -1, num);   %make a spaced vector from 10^-6 to 10^-1

PUE1_log = zeros(size(p_interval));
PUE2_log = zeros(size(p_interval));

% REMEMBER THAT n = 10 (codeword lenght)
for idx= 1:num
    p = p_interval(idx);
    for w = 1:n
    PUE1_log(idx) = PUE1_log(idx) + mult1(w+1) * p^w * (1-p)^(n-w);
    PUE2_log(idx) = PUE2_log(idx) + mult2(w+1) * p^w * (1-p)^(n-w);
    end
end

figure
loglog(p_interval, PUE1_log, 'b-', 'LineWidth', 2)
hold on
loglog(p_interval, PUE2_log, 'r--', 'LineWidth', 2)

grid on
grid minor
xlabel('p')
ylabel('P_{UE}(p)')
title('Undetected-error probability on BSC')
legend('Code 1', 'Code 2', 'Location', 'best')
xlim([1e-6 1e-1])

%% 2.3
% For small p, find the leading term of the two PUE expressions and explain its relationship with d_min and A_min.

%R: For small probabilities, the '(1-p)' part of the sum, can be discarded and the new
%formula contains only the multiplicity part multiplied by probability.
% For linear codes, the dmin between c1 and c2 is equivalent to the result
% of the XOR between the two codewords that is the hamming weight.


%% 2.4
% Compare the two codes: which one has better error-detection performance, and why?

%R: As you can see in the figure the code1 rather than the code2 has a less
%ripid curve that means better performance, for a logarithmic graph the dmin
%leads the gradient because log(Amin) is a constant.
% log Pue(p) = log(Amin) + dmin*log(p)
% for small p (p->0) the formula can be approximate to P_UE ~=
% A_dmin*p^dimn so:
% for G1 the first not null term is for i =4 --> A_min = 15 --> the main
% term P_UE ~= 15p^4
% for G2 the first not null term is for i =3 --> A_min = 1 --> the main
% term P_UE ~= 2p^3



%% part 4

p = [0.10, 0.09, 0.08, 0.07, 0.06];     % given prob for part 4
k = 5;   % information length != codeword length = (info length + redundancy bits)
n = 10; % codeword length

% 4.1 (in the for cycle)


% BSC modeling: y = c_rand + e | where y is the codeword received and "e" are errors during the transmission

P = G1(:, k+1:end);     %parity submatrix (from 6 to 10 columns [5x5])
H = [P', eye(n-k)];    % Syndrome matrix

target_errors = 10;     % stop the cycle at undetected error = 10;
% we need to stop the cycle at 10 undetected errors and then at 100
for j=1:2
    for i=1:length(p)
    
        undetected_errors = 0;  % self-explanatory
        total_transmitted = 0;  % //
        
        while undetected_errors < target_errors
            % 4.1
            u = randi([0, 1], 1, k);     % random generation of a [1 x k] only for 0 and 1
            c_rand = mod(u * G1, 2);    % random gen of a vector 1 x 10 (codeword length)
            
            e = rand(1, length(c_rand)) < p(i);  % generate random noise in BSC with less probability than current
            y = mod(c_rand + e, 2);              % received message
    
            % 4.2 Syndrome calculus
            s = mod(y * H', 2);
    
            % 4.3 Error detection
            if all(s == 0) && ~isequal(y, c_rand)
                undetected_errors = undetected_errors + 1;  % undetected increment
            end
    
            total_transmitted = total_transmitted + 1;  % transmitted anyway
    
        end

    target_errors = 100;    % change the target to 100 errors
    end
end


