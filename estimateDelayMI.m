function [tau, out] = estimateDelayMI(x, maxLag, varargin)
% estimateDelayMI - 用互信息 I[x(t); x(t+lag)] 的第一局部极小估计延迟 tau
% USAGE:
%   [tau, out] = estimateDelayMI(x, maxLag)
%   [tau, out] = estimateDelayMI(x, maxLag, 'Bins', 32, 'Plot', true)
%
% INPUT:
%   x       : 一维列向量
%   maxLag  : 最大滞后（建议 ~ N/10 或按采样物理上限）
%
% NAME-VALUE:
%   'Bins'  : 直方图分箱数（默认 32）
%   'Plot'  : 是否绘图（默认 true）
%
% OUTPUT:
%   tau     : 建议的延迟
%   out     : 结构体，含 I(lag) 序列等

p = inputParser;
p.addParameter('Bins', 32, @isnumeric);
p.addParameter('Plot', true, @islogical);
p.parse(varargin{:});
bins = p.Results.Bins;
doplot = p.Results.Plot;

x = x(:);
x = (x - mean(x)) / std(x);   % 标准化
N = numel(x);
maxLag = min(maxLag, floor(N/2));  % 防止过大

I = zeros(maxLag,1);
edges = linspace(min(x), max(x), bins+1);

for lag = 1:maxLag
    x1 = x(1:end-lag);
    x2 = x(1+lag:end);

    % 单变量直方图
    px = histcounts(x1, edges, 'Normalization','probability');
    py = histcounts(x2, edges, 'Normalization','probability');

    % 联合直方图
    pxy = histcounts2(x1, x2, edges, edges, 'Normalization','probability');

    % 互信息（避免 log(0)）
    px(px==0) = [];
    py(py==0) = [];
    Hx = -sum(px .* log(px));
    Hy = -sum(py .* log(py));

    pxy_vec = pxy(:);
    pxy_nz = pxy_vec(pxy_vec>0);
    Hxy = -sum(pxy_nz .* log(pxy_nz));

    I(lag) = Hx + Hy - Hxy;
end

% 选第一局部极小值；若无，取全局最小（排除 lag=1 的噪声极小可选用平滑）
tau = NaN;
for k = 2:maxLag-1
    if I(k) < I(k-1) && I(k) < I(k+1)
        tau = k;
        break;
    end
end
if isnan(tau)
    [~, tau] = min(I(2:end));  % 避开 lag=1
    tau = tau + 1;
end

out = struct('I', I, 'edges', edges, 'bins', bins);

if doplot
    figure('Color','w');
    plot(1:maxLag, I, 'LineWidth',1.5); grid on;
    xlabel('lag'); ylabel('Mutual Information');
    title(sprintf('Mutual Information vs lag (\\tau = %d)', tau));
    hold on; plot(tau, I(tau), 'o', 'MarkerSize',8);
end
end
