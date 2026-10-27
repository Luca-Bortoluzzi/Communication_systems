clc;
clear;
close all;

%% 1. Code Parameters and Matrix Setup
% Non-systematic generator matrix G (3 x 6)
G = [1 1 0 1 0 1; 
    0 1 1 0 1 1; 
    1 1 1 1 1 0];

k = size(G, 1); % Length of information bits (k = 3)
n = size(G, 2); % Length of codeword bits (n = 6)

% dec2bin(0:2^k-1) generates binary strings for numbers 0 to 7.
% Subtracting '0' (ASCII 48) converts the character array into a numeric binary matrix (8 x 3).
v = dec2bin(0:2^k-1) - '0';  

% Matrix multiplication modulo 2 generates all 2^k possible transmitted codewords (8 x 6).
c = mod(v*G, 2);             

% Generates all 2^n possible received vectors from the channel (64 x 6).
y = dec2bin(0:2^n-1) - '0';  

%% 2. Hamming Distance Matrix Calculation
% Pre-allocate the distance matrix (8 x 64): rows = codewords (c), columns = received vectors (y)
Hd = zeros(2^k, 2^n);

for i = 1:2^k
    for j = 1:2^n
        % xor compares elements bit by bit (returns 1 if different, 0 if equal).
        % sum counts the total number of differing bits (Hamming Distance).
        Hd(i, j) = sum(xor(c(i, :), y(j, :)));
    end
end

%% 3. Minimum Distance Decoding with Random Tie-Breaking
% Array to store the selected codeword index (from 1 to 8) for each received vector y_j

min_val = min(Hd, [], 1);

%% 4. Display Decision Regions D(c_i) with Perfect Alignment
fprintf('=========================================================\n');
fprintf('      DECISION REGIONS WITH RANDOM TIE-BREAKING         \n');
fprintf('=========================================================\n\n');
rule_check=input('select a value from 1 to 8:\n');
flag=NaN;
for i = 1:2^k
    % Logical indexing: select logical 1s where the assigned codeword index equals i
    vector_indices = (Hd(i, :) == min_val);

    % Extract all received vectors y that belong to the decision region D(c_i)
    y_region = y(vector_indices, :);
    if i==1
        D_0=y_region;
    end
    if i==rule_check
        D_C = mod(c(i,:)+D_0,2);
        if sort(y_region)==sort(D_C)
            flag=true;  
        else
            flag=false;
        end
    end

    % sprintf('%d', c(i,:)) formats the binary vector compactly without extra spaces (e.g., [1 0 1 1 1 0] -> 101110)
    c_str = sprintf('%d ', c(i, :)); % Converts vector with a single space between bits

    fprintf('Codeword c(%d) = [%s] --> Contains %d received vectors y:\n', ...
        i, strtrim(c_str), size(y_region, 1));

    % Check if the region contains any vectors before displaying
    if ~isempty(y_region)
        % Loop through each row of y_region and print it with a leading tab (\t)
        for r = 1:size(y_region, 1)
            % Print each bit formatted as %d without extra spaces inserted by num2str
            y_str = sprintf('%d ', y_region(r, :));
            fprintf('\t[%s]\n', strtrim(y_str));
        end
    else
        fprintf('\t(No vectors assigned to this decision region)\n');
    end
    fprintf('---------------------------------------------------------\n');
end
if flag==true
    fprintf("Rule verified\n");
else 
    fprintf('Rule not Verified\n');
end


%% Because of the repetition of codeword in every decision region, there are 10 Codewords, although they're are not the same for every decision region
