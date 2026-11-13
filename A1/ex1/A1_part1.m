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

% Colori coordinati con la presentazione Beamer
navy = [18, 43, 73] / 255;
teal = [0, 128, 128] / 255;

% --- Figura 1: Molteplicità A_i per G1 e G2 (Slide 5) ---
fig1 = figure('Color', 'w', 'Position', [100, 100, 900, 360]);

subplot(1, 2, 1);
b1 = bar(0:10, mult1, 'FaceColor', navy, 'EdgeColor', 'none', 'BarWidth', 0.65);
title('Weight Distribution for G_1', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Hamming weight i', 'FontSize', 11);
ylabel('Multiplicity A_i', 'FontSize', 11);
xticks(0:10);
ylim([0, max(mult1)*1.18]);
grid on; box off;
text(b1.XEndPoints, b1.YEndPoints, string(mult1), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);

subplot(1, 2, 2);
b2 = bar(0:10, mult2, 'FaceColor', teal, 'EdgeColor', 'none', 'BarWidth', 0.65);
title('Weight Distribution for G_2', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Hamming weight i', 'FontSize', 11);
ylabel('Multiplicity A_i', 'FontSize', 11);
xticks(0:10);
ylim([0, max(mult2)*1.18]);
grid on; box off;
text(b2.XEndPoints, b2.YEndPoints, string(mult2), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);

exportgraphics(fig1, 'E1_P1_Q4.png', 'Resolution', 300);

% --- Figura 2: Distanze minime per G1 e G2 (Slide 6) ---
fig2 = figure('Color', 'w', 'Position', [100, 100, 900, 340]);

subplot(1, 2, 1);
h1 = histogram(min1, 'BinMethod', 'integers', 'BinWidth', 0.45, 'FaceColor', navy, 'FaceAlpha', 1);
title('Minimum Distances for G_1', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Minimum Hamming distance', 'FontSize', 11);
ylabel('Frequency', 'FontSize', 11);
xticks(2:5); xlim([1.5, 5.5]);
ylim([0, max(h1.Values)*1.18]);
grid on; box off;
centers1 = h1.BinEdges(1:end-1) + diff(h1.BinEdges)/2;
nz1 = h1.Values > 0;
text(centers1(nz1), h1.Values(nz1), string(h1.Values(nz1)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);

subplot(1, 2, 2);
h2 = histogram(min2, 'BinMethod', 'integers', 'BinWidth', 0.45, 'FaceColor', teal, 'FaceAlpha', 1);
title('Minimum Distances for G_2', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Minimum Hamming distance', 'FontSize', 11);
ylabel('Frequency', 'FontSize', 11);
xticks(2:5); xlim([1.5, 5.5]);
ylim([0, max(h2.Values)*1.18]);
grid on; box off;
centers2 = h2.BinEdges(1:end-1) + diff(h2.BinEdges)/2;
nz2 = h2.Values > 0;
text(centers2(nz2), h2.Values(nz2), string(h2.Values(nz2)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);

exportgraphics(fig2, 'E1_P1_Q5.png', 'Resolution', 300);