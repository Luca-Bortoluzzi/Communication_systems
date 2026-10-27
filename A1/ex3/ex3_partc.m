%% Exercise 3 part C

clc;
clear;
close all;

% To generate a non-linear code, we need to NOT satisfy the property of
% binary sum, two Codewords don't always give as a result another codeword
% of the codebook, and it is non mandatory that the 0_ is present in the
% code

% 3.1 Let's start defining Codewords manually C(4,2)

C = [0 0 0 0 ;
     0 0 1 0 ;
     0 0 1 1 ;
     1 1 0 0];

n = 4; %codeword length
k = 2; %information length

% proving the property for non-linear codes

c1 = C(2,:); % 1 1 1 0 0
c2 = C(4,:); % 0 0 1 1 1
verify = mod(c1 + c2, 2);
fprintf("The codeword sum of the 2nd and the 4th is: %s\n", sprintf("%d", verify));

% this specific xor is not a codeword in our matrix C

% 3.2 Using MATLAB, determine and the decision regions D(c)) for all the four Codewords c ∈ C.

y = dec2bin(0:(2^n-1))- '0'

Hd = zeros(2^k, 2^n);

for i = 1:2^k
    for j = 1:2^n
        % xor compares elements bit by bit (returns 1 if different, 0 if equal).
        % sum counts the total number of differing bits (Hamming Distance).
        Hd(i, j) = sum(xor(C(i, :), y(j, :)));
    end
end

min_val = min(Hd, [], 1);

%% 3.3 Display Decision Regions D(c_i) with Perfect Alignment
fprintf('=========================================================\n');
fprintf('      DECISION REGIONS WITH RANDOM TIE-BREAKING         \n');
fprintf('=========================================================\n\n');

% select a codeword from the given C matrix
while true
    prompt = sprintf("select a value from 1 to %d:\n", n);
    rule_check = input(prompt);
    flag = NaN;
    if(rule_check >= 1 && rule_check <= 4)
        break;
    end
end
for i = 1:2^k
    % Logical indexing: select logical 1s where the assigned codeword index equals i
    vector_indices = (Hd(i, :) == min_val);

    % Extract all received vectors y that belong to the decision region D(c_i)
    y_region = y(vector_indices, :);
    if i == 1
        D_0 = y_region;
    end

    if i == rule_check
        % D(c) = c + D(0)
        D_C = mod(C(i,:) + D_0, 2);
        if isequal(sortrows(y_region), sortrows(D_C))
            flag = true;  
        else
            flag = false;
        end
    end
    % sprintf('%d', c(i,:)) formats the binary vector compactly without extra spaces (e.g., [1 0 1 1 1 0] -> 101110)
    c_str = sprintf('%d ', C(i, :)); % Converts vector with a single space between bits

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

if flag == true
    fprintf("Rule verified\n");
else 
    fprintf('Rule not Verified\n');
end

% 3.4 and 3.5
% The property D(C) = C + D(0) is not verified for non-linear codes, if we
% try inserting the rows of C form the 2nd to the 4th the property is not
% verified, instead for the 1st row because of the fact that the 1st row is
% equal to D(0) for construction, the rule is verified for it is a trivial
% case: D(0) = 0 + D(0)
