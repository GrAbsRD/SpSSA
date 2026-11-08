function [sgc, time] = DRSGMD(x,t,fs,tau,order,threshold_corr,threshold_nmse,sp_components,init_dd,varargin)
% Author: M73ACat
% 2023/03/15
% Basic parameters:
% x                 list, input a signal x,
% fs                int, sample frequency of x,
% tau               int, delay time (optional, defualt 1)
% threshold_corr    float, the threshold of the corr (optional, default 0.95),
% threshold_nmse    float, the threshold of the nmse (optional, default 0.95),
% init_dd           int, manual set the row of the traj-matrix

% Advanced parameters:
% nfft              int, the number of FFT points used to calculate the PSD estimate (optional, default 0 -> 2^floor(log2(n)) ),
% PSD_function      char, functions for calculating PSD (optional ['0'(means 'periodogram'),'1'(means 'pwelch')], default '0' -> periodogram)

% Reference:
% [1] 潘海洋. 基于辛几何模态分解和支持矩阵机的机械故障诊断方法[D]. 湖南大学, 2019.

if nargin < 9
    init_dd = 0;
end
if nargin < 8
    sp_components = [];
end
x = x(:);
n = length(x);

tic;
% trend mode
A = fliplr(vander(t));      % n×n
A = A(:, 1:order+1);            % 取 0~k 次项，n×(k+1)
[U, D, V] = svd(A,'econ');
beta = V * (D \ (U' * x));
trend_component = A * beta;
trend_component = 0.5*t+0.1*t.*t;
trend_component = trend_component(:);
x = x - trend_component;

% symplectic transform
tel=10^-3;
PSD_function = '0';
nfft = 2^floor(log2(n));
for i=1:length(varargin)
    if ischar(varargin{i})
        PSD_function = varargin{i};
    else
        if varargin{i} > 0; nfft = varargin{i}; end
    end
end
if init_dd == 0
    if PSD_function == '0'
        [a,f] = periodogram(x,[],'twosided',nfft,fs);
    else
        window = hanning(nfft);
        [a,f] = pwelch(x,window,nfft*(3/4),nfft,fs,'twoside');
    end
    [~,ra]=max(a(1:end/2));
    if f(ra)/fs<tel
        d=floor(n/3);
    else
        d=floor(1.2*fs/f(ra));
    end
else
    d = init_dd;
end

if ~isempty(sp_components)
    sgc = sp_components;
else
    [X,m]= Trajectorymatrix(x,d,tau);
    A=X'*X;
    [Q, ~] = eig(A);
    Q = real(Q);
    sgc=zeros(n,d);
    for nn=1:d
        z=Q(:,nn)*Q(:,nn)'*X';
        y=zeros(n,1);
        newd=min(size(z));
        newm=max(size(z));
        if m>d
            z1=rot90(z',3);
        else
            z1=rot90(z,3);
        end
        for k=1:n
            if k<=newd
                y(k)=sum(diag(z1,newm-k))/k;
            elseif k>newd&&k<=newm
                y(k)=sum(diag(z1,newm-k))/newd;
            else
                y(k)=sum(diag(z1,newm-k))/(n-k+1);
            end
        end
        sgc(:,nn)=y;
    end
end
time = toc;
if std(sgc(:,end)) > std(sgc(:,1))
    y = fliplr(sgc);
else
    y = sgc;
end

% noise mode
if threshold_nmse == 0
    noise = zeros(size(x));
    d_ = d;
else
    result = noiseColsByColNormEnergy(y, 1-threshold_nmse);
    noise = sum(y(:, result.noise_idx), 2);
    y = y(:, result.keep_idx);
    d_ = length(result.keep_idx);
end

% oscillatory mode
sgc = trend_component;
alive = true(1, d_);     % 还没被分组的列
while any(alive)
    i = find(alive, 1, 'first');   % 取第一个还活着的索引
    source = y(:, i);
    alive(i) = false;
    cand = find(alive);            % 当前剩余候选
    if ~isempty(cand)
        keep = false(size(cand));  % 记录要合并进来者
        for k = 1:numel(cand)
            j = cand(k);
            sim1 = simFreqDiffByHilbertFFT(source, y(:, j), fs);
            sim2 = cosSim(source, y(:, j));
            % sim2 = TGSI(source, y(:, j));
            if sim1 > threshold_corr || sim2 > threshold_corr
                keep(k) = true;
            end
        end
        grp = cand(keep);          % 这一轮与 source 相似的列
        if ~isempty(grp)
            alive(grp) = false;
            source = source + sum(y(:, grp), 2);  % 直接累加
        end
    end
    sgc(:,end+1) = source; %#ok<AGROW>
end
sgc(:,end+1) = noise;
end

function [t,m] = Trajectorymatrix(x,d,tau)
%根据Takens嵌入定理构建的轨迹矩阵
%input
%x          input a signal x,
%d          embedding dimension
%tau        delay time
%output
%t m*n Trajectory matrix
%20180611 PhD. Jin Hang   in Chengdu
[n ,c]=size(x);
if n==1 || c==1
    if c>n
        x=x';
        [n ,~]=size(x);
    end
else
    error('x isn`t  Column vector')
end
m=n-(d-1)*tau;
if m<=0
    error('Index out of matrix,Replace d or Tau')
end
t=zeros(m,d);
for kk=1:d
    t(1:m,kk)=x(1+(kk-1)*tau:m+(kk-1)*tau);
end
end

function result = noiseColsByColNormEnergy(X, keepRatio)
% removeNoiseAndSortCols
% 通过列能量 E_j = ||x_j||_2^2 选取有效列，剩余视为噪声列；
% 并将去噪后的矩阵按能量从大到小重排。
%
% 输入:
%   X         : n×m 矩阵（每列一个向量）
%   keepRatio : (0,1] 保留能量比例，比如 0.95
%
% 输出:
%   noise_idx      : 被判为噪声的列索引（相对原矩阵）
%   keep_idx_sorted: 保留下来的列索引（相对原矩阵，且按能量降序）

    if nargin < 2 || isempty(keepRatio), keepRatio = 0.95; end

    E = sum(X.^2, 1);                      % 1×m, 每列能量
    [~, order] = sort(E, 'descend');       % 能量降序的原列索引
    totalE = sum(E);

    cumRatio = cumsum(E(order)) / totalE;  % 降序累计能量占比
    r = find(cumRatio >= keepRatio, 1, 'first');
    if isempty(r), r = numel(order); end

    keep_idx_sorted = order(1:r);
    noise_idx = sort(order(r+1:end));      % 噪声索引升序便于查看

    result.noise_idx = noise_idx;
    result.keep_idx = keep_idx_sorted;

end

function [omega0k, f0k] = centerFreq_argmax(y, fs)
% centerFreq_argmax: ω0k = argmax |F{ y + j*H{y} }|
% y        : 实信号 (N×1 或 1×N)
% fs       : 采样率 Hz
% useWindow: 是否加窗抑制泄漏（可选，默认 true）
%
% omega0k  : 角频率 rad/s
% f0k      : 频率 Hz
% fAxis    : 单边频率轴 Hz
% magSpec  : 单边幅度谱 |FFT(analytic)|

    y = y(:) - mean(y);
    % Hilbert 包络
    env = abs(hilbert(y));
    
    % FFT
    N = length(env);
    Y = abs(fft(env));
    f = (50:N-1) * fs / N;
    
    % 4) use one-sided spectrum (analytic signal ideally has +freq only)
    half = 1:floor(N/2);

    % 5) argmax
    [~, idx] = max(Y(half));
    f0k = f(idx);
    omega0k = 2*pi*f0k;
end

function [sim1, omega1, omegak] = simFreqDiffByHilbertFFT(y1, yk, fs)

    omega1 = estimate_main_peak_3gauss(y1, fs);
    omegak = estimate_main_peak_3gauss(yk, fs);

    % --- 4) 频率差比 sim1 ---
    sim1 = 1 - abs(omega1 - omegak) / max(omega1, omegak);

end

function c = cosSim(y1, y2)
    y1 = y1(:); y2 = y2(:);
    N = min(numel(y1), numel(y2));
    y1 = y1(1:N); y2 = y2(1:N);

    den = norm(y1) * norm(y2);
    if den == 0
        c = NaN; % 或者 0
    else
        c = real((y1' * y2) / den);  % 复信号取 real，避免数值虚部
    end
end
