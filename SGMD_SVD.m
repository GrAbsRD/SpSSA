
% SymplecticSingularSpecAnalysis
function [SPSSCs, time] = SGMD_SVD( ...
    sig, fs, tau, n_cluster, ...
    threshold_corr, threshold_nmse, ...
    merge_mode, periodic_extension, init_dd ...
    )

if nargin < 9
    init_dd = 0;
end
if nargin < 8
    periodic_extension = false;
end
if nargin < 7
    merge_mode = 'kmeans';
end

tic
%% embedding
n = length(sig);

if init_dd == 0
    tel = 1e-3;
    nfft = 2^floor(log2(n));
    [psd,freqs] = periodogram(sig, [], 'twosided', nfft, fs);
    [~, ra] = max(psd(1:floor(end/2)));
    if freqs(ra) / fs < tel
        d = floor(n / 3);
    else
        d = floor(1.2 * fs / freqs(ra));
    end
else
    d = init_dd;
end
m = n - (d - 1) * tau;
if periodic_extension
    offset = (d - 1) * tau;
    f = @(k) mod(k-1, n) + 1;
    sig = sig(f(1:n+offset));
end
trajectory_matrix = TrajectoryMatrix(sig, d, tau);

if d >= m
    error("d must less than m !")
end

A = trajectory_matrix * trajectory_matrix';
[Q, ~] = eig(A * A);
% [Q, ~] = eig(A);
Q = real(Q);
% [U, D, V] = svd(trajectory_matrix);

%% pull back
symplectic_component = zeros(n, d);
for i = 1:d
    reconstructed_trajectory_matrix = Q(:, i)*Q(:, i)'*trajectory_matrix;
    % reconstructed_trajectory_matrix = D(i, i) * U(:, i) * V(:, i)';
    new_d = min(m, d);
    new_m = max(m, d);
    if d > m
        reconstructed_trajectory_matrix = reconstructed_trajectory_matrix';
    end
    if periodic_extension
        reconstructed_trajectory_matrix = [
            reconstructed_trajectory_matrix(:, m+1:n), ...
            reconstructed_trajectory_matrix ...
            ];
        nn = n + offset;
        new_m = new_m + 2 * offset;  % size(reconstructed_trajectory_matrix, 2)
    else
        nn = n;
    end
    % 生成网格索引
    [P, J] = ndgrid(1:new_d, 1:new_m);          % P: 行索引 p, J: 列索引 j
    K = J + P*tau - tau;                        % k = j + p*tau - tau
    mask = (K >= 1) & (K <= nn);
    % 取出所有有效元素
    vals = reconstructed_trajectory_matrix(mask);
    kk   = K(mask);
    % 对每个 k 做 sum 和 count
    s = accumarray(kk(:), vals(:), [nn 1], @sum, 0);
    c = accumarray(kk(:), 1,       [nn 1], @sum, 0); 
    single_component = s ./ c;
    single_component(c==0) = 0;   % 若有空（理论上不应太多），避免 NaN
    % single_component = zeros(nn, 1);
    % if tau == 1 && 0
    %     z = rot90(reconstructed_trajectory_matrix', 3);
    %     for k = 1:nn
    %         single_component(k) = mean(diag(z, new_m-k));
    %     end
    % else
    %     for k = 1:nn
    %         left = max(1, ceil((k + tau - new_m) / tau));
    %         right = min(new_d, floor((k + tau - 1) / tau));
    %         tmp = zeros(right - left + 1, 1);
    %         for p = left:right
    %             tmp(p - left + 1) = ...
    %                 reconstructed_trajectory_matrix(p, k-p*tau+tau);
    %         end
    %         single_component(k) = mean(tmp);
    %     end
    % end
    if periodic_extension
        single_component = single_component(offset+1:end);
    end
    symplectic_component(:, i) = single_component;
end
time = toc;

if std(symplectic_component(:,end)) > std(symplectic_component(:,1))
    symplectic_component = fliplr(symplectic_component);
end

%% merging
if strcmp(merge_mode, 'kmeans')
    if n_cluster == 1
        reconstructed_signal = sum(symplectic_component, 2);
    else
        [labels, ~] = kmeans(symplectic_component', n_cluster, ...
            'Replicates', 5);
        tmp = unique(labels, 'stable');
        reconstructed_signal = zeros(n, n_cluster);
        for i = 1:n_cluster
            single_signal = sum(symplectic_component(:, labels==tmp(i)), 2);
            reconstructed_signal(:, i) = single_signal(1:n);
        end
    end
    existing_signal = sum(symplectic_component, 2);
    % added noise
    reconstructed_signal = [reconstructed_signal, sig(1:n) - existing_signal(1:n)];
else
    cosSim = @(x,y) (x(:)' * y(:)) / (norm(x) * norm(y));
    reconstructed_signal = [];
    indexes = 1:d;
    if strcmp(merge_mode, 'std')
        x_e = std(sig);
    elseif strcmp(merge_mode, 'var')
        x_e = var(sig);
    else  % nmse
        x_e = sum((sig - mean(sig)).^2);
    end
    remaining_signal = sig;
    while ~isempty(indexes)
        flags = logical(indexes);
        source = symplectic_component(:, indexes(1));
        flags(1) = 0;
        grouped_idx = [indexes(1)];
        for i = 2:length(indexes)
            idx = indexes(i);
            corrs = corrcoef(source, symplectic_component(:, idx)); 
            corrs = corrs(1, 2);
            % corrs = norm(source - symplectic_component(:, idx));
            if corrs >= threshold_corr
                grouped_idx(end + 1) = idx;
                flags(i) = 0;
            end
        end
        indexes = indexes(flags == 1);
        single_signal = sum(symplectic_component(:, grouped_idx), 2);
        reconstructed_signal = [reconstructed_signal, single_signal];
        existing_signal = sum(reconstructed_signal, 2);
        remaining_signal = sig(1:n) - existing_signal(1:n);
        if strcmp(merge_mode, 'std')
            g_e = std(remaining_signal);
        elseif strcmp(merge_mode, 'var')
            g_e = var(remaining_signal);
        else  % nmse
            g_e = sum(remaining_signal.^2);
        end
        if g_e / x_e < threshold_nmse
            break
        end
    end
    reconstructed_signal = [reconstructed_signal, remaining_signal];  % noise
end

SPSSCs = reconstructed_signal;

end


%% Auxiliary function
function M = TrajectoryMatrix(x, row, tau)
n = length(x);
m = n - (row - 1) * tau;
M = zeros(row, m);
for k = 1:row
    offset = (k - 1) * tau;
    M(k, 1:m) = x(1+offset:m+offset);
end
%M = [M(:, 1+n-2*(row-1)*tau:m) M(:, 1:n-2*(row-1)*tau)];
end

function XY = DiophantineEqSolver(n, tau, x1, x2, y1, y2, offset_x)
left = max(x1, ceil((n - tau * offset_x - y2) / tau));
right = min(x2, floor((n - tau * offset_x - y1) / tau));
XY = zeros(right-left+1, 2);
k = 1;
for x = left:right
    XY(k, :) = [x, n - tau * (x + offset_x)];
    k = k + 1;
end
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
