function [m_opt, out] = estimateEmbeddingFNN(x, tau, mMax, varargin)
% estimateEmbeddingFNN - 用 False Nearest Neighbors 估计嵌入维数 m
% 参考：Kennel et al. (1992)
%
% USAGE:
%   [m_opt, out] = estimateEmbeddingFNN(x, tau, mMax)
%   [m_opt, out] = estimateEmbeddingFNN(x, tau, mMax, ...
%                   'Rtol', 15, 'Atol', 2, 'W', tau, 'Plot', true, 'Thresh', 0.01)
%
% INPUT:
%   x     : 一维列向量
%   tau   : 延迟
%   mMax  : 最大嵌入维
%
% NAME-VALUE:
%   'Rtol'  : 距离比阈值（默认 15）
%   'Atol'  : 吸引子尺度阈值（单位：std(x)，默认 2）
%   'W'     : Theiler window（时间上剔除近邻，默认 tau）
%   'Plot'  : 是否绘图（默认 true）
%   'Thresh': FNN 比例阈值（默认 0.01，即 1%）
%
% OUTPUT:
%   m_opt : 建议的嵌入维（FNN% 低于阈值的最小 m；若无则给出肘点）
%   out   : 结构体，含 FNN% 曲线等

p = inputParser;
p.addParameter('Rtol', 15, @isnumeric);
p.addParameter('Atol', 2, @isnumeric);
p.addParameter('W', [], @isnumeric);
p.addParameter('Plot', true, @islogical);
p.addParameter('Thresh', 0.01, @isnumeric);
p.parse(varargin{:});
Rtol = p.Results.Rtol;
Atol = p.Results.Atol;
W = p.Results.W;
doplot = p.Results.Plot;
thresh = p.Results.Thresh;

x = x(:);
x = (x - mean(x)) / std(x);
N = numel(x);
if isempty(W), W = tau; end

% 构造嵌入向量函数
embed = @(m) hankel_embed(x, m, tau);   % 返回 M×m

fnn_ratio = nan(mMax,1);
fnn_count = nan(mMax,1);
total_count = nan(mMax,1);

for m = 1:mMax
    X_m = embed(m);
    M = size(X_m,1);
    if M < 2
        fnn_ratio(m) = NaN; continue;
    end

    % 最近邻（排除 Theiler window）
    fnn = 0; tot = 0;

    % 为提速，使用 knnsearch
    IDX = knnsearch(X_m, X_m, 'K', 2);  % 自身+最近邻
    nn  = IDX(:,2);                     % 最近邻索引

    for i = 1:M
        if abs(i - nn(i)) <= W
            % 近时间的邻居弃用，找下一个可用邻居
            % 简易策略：扩大K再找一个可用邻居
            k = 3;
            while true
                if k > min(M,10) % 最多尝试到10近邻
                    nn_i = NaN; break;
                end
                idxk = knnsearch(X_m(i,:), X_m, 'K', k);
                idxk = idxk(2:end); % 去掉自身
                idxk = idxk(abs(idxk - i) > W);
                if ~isempty(idxk)
                    nn_i = idxk(1); break;
                end
                k = k + 1;
            end
            if isnan(nn_i), continue; end
        else
            nn_i = nn(i);
        end

        % m 维距离
        dist_m = norm(X_m(i,:) - X_m(nn_i,:));
        if ~isfinite(dist_m) || dist_m==0
            continue;
        end

        % m+1 维距离（需要 m+1 维嵌入）
        if m < mMax
            X_m1 = embed(m+1);
            if i > size(X_m1,1) || nn_i > size(X_m1,1), continue; end
            dist_m1 = norm(X_m1(i,:) - X_m1(nn_i,:));
        else
            % 最后一维就不计算 m+1 的判据了
            continue;
        end

        % FNN 条件（Kennel 两个判据）
        cond1 = sqrt(max(dist_m1^2 - dist_m^2, 0)) / dist_m > Rtol;
        cond2 = abs( X_m1(i,end) - X_m1(nn_i,end) ) > Atol; % 因为 x 已 z-score

        if cond1 || cond2
            fnn = fnn + 1;
        end
        tot = tot + 1;
    end

    fnn_count(m) = fnn;
    total_count(m) = tot;
    fnn_ratio(m) = (tot>0) * (fnn / max(1,tot));
end

% 选择 m_opt
m_opt = find(fnn_ratio <= thresh, 1, 'first');
if isempty(m_opt)
    % 若没有低于阈值，取“肘点”：差分最小处
    dff = diff(fnn_ratio);
    [~, idx] = min(dff);
    m_opt = idx + 1;
end

out = struct('fnn_ratio', fnn_ratio, 'fnn_count', fnn_count, ...
             'total', total_count, 'Rtol', Rtol, 'Atol', Atol, ...
             'W', W, 'thresh', thresh);

if doplot
    figure('Color','w');
    stem(1:mMax, fnn_ratio, 'filled'); grid on;
    xlabel('Embedding dimension m'); ylabel('FNN proportion');
    title(sprintf('FNN vs m (m_{opt}=%d, thresh=%.2g)', m_opt, thresh));
    hold on; plot(m_opt, fnn_ratio(m_opt), 'o', 'MarkerSize',8);
end
end

% ------- helper: 构造延迟嵌入矩阵 -------
function Xe = hankel_embed(x, m, tau)
N = numel(x);
M = N - (m-1)*tau;
if M <= 0, Xe = zeros(0,m); return; end
Xe = zeros(M, m);
for j = 1:m
    Xe(:,j) = x( (1:M) + (j-1)*tau );
end
end
