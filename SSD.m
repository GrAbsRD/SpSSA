function [components, resid, info] = SSD(sig, fs, tau, opts)
%SSD_REFACTORED Singular Spectrum Decomposition (Bonizzi et al. 2014)
%
% INPUTS:
%   sig  - Input signal (vector)
%   fs   - Sampling frequency (Hz)
%   opts - (Optional) Struct with fields:
%          .threshold_nmse : Stop if NMSE < thr (default 0.01 [cite: 230])
%          .max_iter       : Maximum iterations (default 50)
%          .trend_flimit   : Normalized freq limit for trend (default 1e-3 [cite: 175])
%
% OUTPUTS:
%   components - Extracted components [N x k]
%   resid      - Final residual signal
%   info       - Diagnostic info struct

    %% 1. Parameter Initialization
    if nargin < 4, opts = struct(); end
    if ~isfield(opts, 'threshold_nmse'), opts.threshold_nmse = 0.01; end
    if ~isfield(opts, 'max_iter'), opts.max_iter = 50; end
    if ~isfield(opts, 'trend_flimit'), opts.trend_flimit = 1e-3; end

    sig = sig(:);
    N = length(sig);
    v = sig; % Residual initialization [cite: 173]
    
    components = [];
    energy_orig = sum(sig.^2);
    
    info.iterations = [];

    %% 2. Iterative Decomposition
    for iter = 1:opts.max_iter
        
        % Check stopping criterion (NMSE) [cite: 229-231]
        current_nmse = sum(v.^2) / energy_orig;
        if current_nmse <= opts.threshold_nmse
            break;
        end
        
        %% 2.1 Estimate Dominant Frequency & Bandwidth [cite: 173-174]
        % (Assuming estimate_main_peak_3gauss is available from context)
        [fmax, delta_f] = estimate_main_peak_3gauss(v, fs);
        
        %% 2.2 Adaptive Window Length M [cite: 175-178]
        norm_f = fmax / fs;
        is_trend = false;
        
        if iter == 1 && norm_f < opts.trend_flimit
            % Trend detection at first iteration
            d = floor(N / 3); 
            is_trend = true;
        else
            % Oscillatory component
            d = floor(1.2 * fs / fmax);
        end
        
        % Safety check for M
        d = max(2, min(floor(N/2), d)); 
        
        %% 2.3 Build Trajectory Matrix (Wrapped) [cite: 147-154]
        % Note: If is_trend is true, Bonizzi suggests standard matrix [cite: 165-168],
        % otherwise use the new wrapped matrix.
        % if is_trend
        %     X = BuildTrajectoryMatrix_Standard(v, d);
        % else
        %     X = BuildTrajectoryMatrix_Wrapped(v, d);
        % end
        if is_trend
            offset = 0;
        else
            offset = (d - 1) * tau;
        end
        X = TrajectoryMatrix(v, d, tau, offset);
        
        %% 2.4 SVD Decomposition [cite: 105-106]
        % Using 'econ' ensures we don't compute unnecessary V columns if M > N
        [U, S, V] = svd(X, 'econ'); 
        s_vals = diag(S);
        d_ = EfficientRank(s_vals, .25);
        
        %% 2.5 Component Selection [cite: 184-189]
        if is_trend
            % For sizable trend, use only the first eigentriple [cite: 184]
            indices = 1;
        else
            % Identify components in the frequency band [fmax-df, fmax+df]
            indices = SelectOscillatoryComponents(U(:, 1:d_), fs, fmax, delta_f);
            % indices = group_components_by_left_sv(U(:, 1:d_), fs, fmax, delta_f);
        end
        
        if isempty(indices)
            % If no component found in band, force select the one contributing 
            % most to the dominant peak (Logic implied in [cite: 188])
            indices = 1; 
        end
        
        %% 2.6 Reconstruction (Diagonal Averaging) [cite: 189]
        % Reconstruct matrix subset: X_I = sum(s_i * u_i * v_i')
        X_sub = U(:, indices) * S(indices, indices) * V(:, indices)';
        
        % if is_trend
        %     g_est = DiagonalAveraging_Standard(X_sub);
        % else
        %     g_est = DiagonalAveraging_Wrapped(X_sub, N); 
        % end
        g_est = GeneralizedDiagonalAveraging(X_sub, N, tau, offset);
        
        %% 2.7 Optimization (Scaling) [cite: 223-224]
        % a_hat = (g^T * v) / (g^T * g)
        num = g_est' * v;
        den = g_est' * g_est;
        if den < 1e-10, den = 1e-10; end % numerical stability
        a_hat = num / den;
        
        g_final = a_hat * g_est;
        % g_final = g_est;
        
        %% 2.8 Update Residual [cite: 228]
        v = v - g_final;
        
        % Store result
        components = [components, g_final];
        
        % Save info
        iter_info.fmax = fmax;
        iter_info.d = d;
        iter_info.is_trend = is_trend;
        iter_info.NMSE = current_nmse;
        info.iterations = [info.iterations; iter_info];
    end
    
    resid = v;
    components = [components, resid];

end

%% ---------------- Helper Functions ---------------- %%

function X = BuildTrajectoryMatrix_Wrapped(x, M)
    % Implements the "wrapped" trajectory matrix from 
    % Rows are cyclic shifts of the signal.
    N = length(x);
    X = zeros(M, N);
    for i = 1:M
        % Row i: x(i), ..., x(N), x(1), ..., x(i-1)
        % This is equivalent to a circular shift
        X(i, :) = circshift(x, -(i-1))'; 
    end
end

function X = BuildTrajectoryMatrix_Standard(x, M)
    % Standard SSA trajectory matrix (Hankel) 
    N = length(x);
    K = N - M + 1;
    X = zeros(M, K);
    for i = 1:M
        X(i, :) = x(i : i+K-1)';
    end
end

function M = TrajectoryMatrix(x, row, tau, offset)  
    if nargin < 4
        offset = 0;
    end
    n = length(x);
    f = @(k) mod(k-1, n) + 1;
    x = x(f(1:n+offset));
    n = n + offset;
    m = n - (row - 1) * tau;
    M = zeros(row, m);
    for k = 1:row
        s = (k - 1) * tau;
        M(k, 1:m) = x(1+s:m+s);
    end
end

function g = GeneralizedDiagonalAveraging(X, N, tau, offset)
    if nargin < 4
        offset = 0;
    end
    if offset > 0
        X = [X(:, end-offset+1:end), X];
    end
    [d, m] = size(X);
    g = zeros(N+offset, 1);
    for k = 1:N+offset
        left = max(1, ceil((k + tau - m) / tau));
        right = min(d, floor((k + tau - 1) / tau));
        tmp = zeros(right - left + 1, 1);
        for p = left:right
            tmp(p - left + 1) = X(p, k-p*tau+tau);
        end
        g(k) = mean(tmp);
    end
    g = g(offset+1:end);
end

function g = DiagonalAveraging_Wrapped(X, N)
    % Adapts diagonal averaging for the wrapped matrix 
    % For the wrapped matrix definition X(i,j) = x((i+j-2 mod N) + 1)
    % The average is along the "extended" cross-diagonals.
    % Actually, for the specific structure in, row i is shifted by i-1.
    % So column j of row i corresponds to time index: mod(j + i - 2, N) + 1.
    % The simplest reconstruction is averaging the "un-shifted" rows.
    
    [M, ~] = size(X);
    g_accum = zeros(N, 1);
    
    for i = 1:M
        % Inverse of the construction shift
        row_i = X(i, :);
        row_i_unshifted = circshift(row_i, (i-1)); 
        g_accum = g_accum + row_i_unshifted(:);
    end
    
    g = g_accum / M;
end

function g = DiagonalAveraging_Standard(X)
    % Standard Hankel Diagonal Averaging [cite: 114]
    [M, K] = size(X);
    N = M + K - 1;
    g = zeros(N, 1);
    
    % Simple implementation of diagonal averaging
    L = min(M, K);
    for k = 1:N
        s = 0;
        count = 0;
        % Determine range of i
        start_i = max(1, k - K + 1);
        end_i = min(k, M);
        
        for i = start_i:end_i
            j = k - i + 1;
            s = s + X(i, j);
            count = count + 1;
        end
        g(k) = s / count;
    end
end

function idx = SelectOscillatoryComponents(U, fs, fmax, delta_f)
    % Select columns of U with dominant frequency in band [cite: 188]
    [M, r] = size(U);
    idx = [];
    
    for i = 1:r
        u_vec = U(:, i);
        u_vec = u_vec - mean(u_vec);
        
        % Simple Peak estimate (can use periodogram or pwelch)
        [Pxx, F] = periodogram(u_vec, [], [], fs);
        [~, max_loc] = max(Pxx);
        dom_freq = F(max_loc);
        
        if dom_freq >= (fmax - delta_f) && dom_freq <= (fmax + delta_f)
            idx = [idx, i];
        end
    end
    
    % Paper implies selecting at least the most energetic one if multiple exist,
    % or pairing them. Here we return all in band.
end

function r = EfficientRank(singular_vec, gamma)
if nargin < 2
    gamma = .97;
end
r = length(singular_vec);
power = 0;
total_power = norm(singular_vec, 2)^2;
for i = r:-1:1
    power = power + singular_vec(i)^2;
    if (total_power - power) / total_power < gamma
        r = i;
        break
    end
end
end