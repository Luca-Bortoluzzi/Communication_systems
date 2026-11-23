clc;
clear;
close all;

%% 18. Generate the 128 x 128 Walsh-Hadamard matrix
N = 128;
K = 128;
Had = hadamard(N);

%% 19. Verify zero-delay orthogonality between distinct rows
is_orth = true;
for i = 1:K-1
    for j = i+1:K
        if (Had(i, :) * Had(j, :)') ~= 0
            is_orth = false;
        end
    end
end

if is_orth
    fprintf("19. Zero-delay orthogonality between distinct rows verified.\n");
else
    fprintf("19. Orthogonality NOT verified.\n");
end

%% 20. Compute r_A, r_C and r_M using cyclic correlations
H_fft = fft(Had, [], 2);
r_A = 0;
r_C = 0;

for i = 1:K
    for j = i:K
        R_periodic = round(real(ifft(H_fft(i, :) .* conj(H_fft(j, :)))));
        
        if i == j
            R_periodic(1) = 0; % Exclude zero-delay peak (N = 128)
            r_A = max(r_A, max(abs(R_periodic)));
        else
            r_C = max(r_C, max(abs(R_periodic)));
        end
    end
end

r_M = max(r_A, r_C);
fprintf("20. Cyclic correlation results: r_A = %d, r_C = %d, r_M = %d\n", r_A, r_C, r_M);

%% 21. Compare with Welch and Sidelnikov bounds for N = K = 128
s_welch = 1:10;
temp_w = zeros(size(s_welch));
for idx = 1:length(s_welch)
    s = s_welch(idx);
    binom = nchoosek(N + s - 1, s);
    arg = max(0, (1 / (K*N - 1)) * ((K*N) / binom - 1));
    temp_w(idx) = N * (arg ^ (1 / (2*s)));
end
welch_gen    = ceil(max(temp_w));
welch_s1     = ceil(N * sqrt((K - 1) / (K*N - 1)));
welch_approx = ceil(sqrt(N));

s_sidelnikov = 0:min(10, ceil(2*N/5) - 1);
temp_sid = zeros(size(s_sidelnikov));
for idx = 1:length(s_sidelnikov)
    s = s_sidelnikov(idx);
    term1 = (2*s + 1) * (N - s);
    term2 = (s * (s + 1)) / 2;
    term3 = (2^s * N^(2*s + 1)) / (K * factorial(2*s) * nchoosek(N, s));
    inner_arg = max(0, ceil(term1 + term2 - term3));
    temp_sid(idx) = ceil(sqrt(inner_arg));
end
sidelnikov = max(temp_sid);

fprintf("21. Bounds at (N=128, K=128): Welch gen = %d | Welch (s=1) = %d | Welch approx = %d | Sidelnikov = %d\n", ...
    welch_gen, welch_s1, welch_approx, sidelnikov);
fprintf("    Hadamard r_M (%d) is far worse than Sidelnikov (%d) and Welch (%d) bounds.\n", ...
    r_M, sidelnikov, welch_gen);

%% 22. Comment on the results, Downlink/Uplink use cases, and Synchronization issues
%
% - Results comparison:
%   For the 128x128 Walsh-Hadamard family, we obtain r_A = 128, r_C = 128, and
%   thus r_M = 128. Although distinct sequences have zero cross-correlation at
%   zero delay (tau = 0), their cyclic correlation for non-zero delays reaches
%   the maximum possible value N = 128 (e.g., the all-ones first row has an
%   autocorrelation of 128 for every cyclic shift, and shifting certain rows
%   aligns them with other rows). Consequently, r_M = 128 is far above the
%   Sidelnikov (16) and Welch (12) lower bounds, making Hadamard codes very poor
%   in terms of periodic/asynchronous correlation.
%
% - Downlink vs. Uplink use cases:
%   * Downlink (Base Station to Mobile Users): Walsh-Hadamard sequences are
%     ideal for synchronous downlink channels (such as WCDMA/CDMA2000). Since all
%     user signals originate from a single transmitter (the base station), they
%     can be transmitted simultaneously with zero relative delay (tau = 0),
%     preserving strict orthogonality and eliminating Multi-User Interference (MUI).
%   * Uplink (Mobile Users to Base Station): They are unsuitable on their own for
%     the uplink, because users are geographically dispersed at different
%     distances from the base station and have independent local clocks, causing
%     signals to arrive asynchronously (tau != 0).
%
% - Synchronization issues:
%   Even small timing offsets, clock jitter, or multipath propagation delays
%   introduce non-zero cyclic/aperiodic shifts (tau != 0). Because r_A = r_C = 128,
%   any loss of synchronization destroys orthogonality and produces severe
%   self-interference and cross-user interference. In practical cellular systems,
%   this is mitigated by multiplying Walsh-Hadamard channelization codes with a
%   long pseudo-noise (PN) scrambling sequence to improve out-of-phase correlation.