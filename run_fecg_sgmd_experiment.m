%% 0) 主脚本
function sgmd_components = run_fecg_sgmd_experiment(x)
clc; close all;

if nargin < 1
    x = [];
end

% ========= 用户参数 =========
csvFile = 'aFECGdata/r01.csv';
Fs = 1000;
useNotch50 = true;
notchHz = 50;

% ========= 1) 读取数据 =========
D = read_fecg_csv(csvFile);
abd = D.Abdomen_1;

% ========= 2) 预处理 =========
% abd_f = preprocess_ecg(abd, Fs, useNotch50, notchHz);
abd_f = build_weighted_sum_GSVD(D, Fs);
% abd_f = (-3*D.Abdomen_1+D.Abdomen_2+D.Abdomen_3+D.Abdomen_4)/3;
% abd_f = abd;

% ========= 3) SGMD 分解 =========
tau = 1;
n_cluster = 2;
threshold_corr = 0.85;
threshold_nmse = 0.001;
threshold_svd = 0.999;
periodic_extension = false;
% sgmd_components: 每一列是一个分解分量
if ~isempty(x)
    sgmd_components = x;
else
    sgmd_components = cell(7, 1);
    % opts.threshold_nmse=0.01;
    % sgmd_components{1} = SSD(abd_f(:),Fs,3,opts);
    % sgmd_components{1} = JMD(abd_f(:), )
    sgmd_components{1} = SGMD(abd_f(:)', Fs, tau, threshold_corr, threshold_nmse);
    sgmd_components{1} = vmd(abd_f(:));
    for tau = 1:3
        sgmd_components{2+(tau-1)*2} = SpSSA(abd_f(:), Fs, tau, n_cluster, ...
            threshold_corr, threshold_nmse, threshold_svd, ...
            'kmeans', false);
        sgmd_components{3+(tau-1)*2} = SpSSA(abd_f(:), Fs, tau, n_cluster, ...
            threshold_corr, threshold_nmse, threshold_svd, ...
            'kmeans', true);
    end
end
% sgmd_components{7} = SpSSA(abd_f(:), Fs, 3, n_cluster, ...
%             threshold_corr, threshold_nmse, threshold_svd, ...
%             'kmeans', true);

%%
fs = 1000;          % 采样率
Tmax = 10;          % 前 10 秒
% 读取表格
T = readtable('aFECGdata/r01_annotations.csv');
% 提取时间和标注
onset = T.Onset;            % 秒
label = T.Annotation;       % 标注字符串
% 1️⃣ 只保留 QRS（不含 R-peak@@Direct_1）
idx_qrs = strcmp(label, 'QRS');
onset_qrs = onset(idx_qrs);
% 2️⃣ 只取前 10 秒
idx_10s = onset_qrs <= Tmax;
onset_qrs_10s = onset_qrs(idx_10s);
% 3️⃣ 转换为样本位置（从 1 开始）
qrs_sample = round(onset_qrs_10s * fs) + 1;

%%
t = linspace(0,10,Fs*Tmax);
fig=figure('Units','pixels', 'Position', [100, 100, 1200, 550]);
tlo = tiledlayout(fig,2,1, 'TileSpacing','compact', 'Padding','compact');
ax = nexttile(tlo, 1);
xline(onset_qrs_10s, '--', 'LineWidth', .5,'Color','red');hold on;
plot(t,abd_f,'LineWidth', .5, 'Color','blue');
ylabel('Amplitude(μV)', 'Interpreter', 'latex','FontSize',18)
p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
p(1) = -0.45;                % 统一的 x 偏移（越小越靠左）
p(2) = 0;                  % 垂直居中（0~1）
ax.YLabel.Position = p;
set(ax,'fontsize',16,'LineWidth',1.0,'TickLength',[0.01, 0.02],'Box','on')
ax = nexttile(tlo, 2);
xline(onset_qrs_10s, '--', 'LineWidth', .5,'Color','red');hold on;
plot(t,D.Direct_1,'LineWidth', .5, 'Color','blue');
xlabel('time(s)', 'FontSize',18)
ylabel('Amplitude(μV)', 'Interpreter', 'latex','FontSize',18)
p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
p(1) = -0.45;                % 统一的 x 偏移（越小越靠左）
p(2) = 0;                  % 垂直居中（0~1）
ax.YLabel.Position = p;
set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02],'Box','on')
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
% print('./pic/fecg_signals.png', '-dpng', '-r600');

%%
% for ii = 1:7
%     fig = figure('Units','pixels', 'Position', [100, 100, 1200, 1050]);
%     tlo = tiledlayout(fig,8,2, 'TileSpacing','compact', 'Padding','compact');  % 关键
%     for i = 1:8
%         for j = 1:2
%             if (i-1)*2+j > size(sgmd_components{ii}, 2)
%                 break
%             end
%             x = sgmd_components{ii}(:, (i-1)*2+j);
%             nexttile(tlo, (i-1)*2+j);
%             % xline(onset_qrs_10s, 'r--', 'LineWidth', 1.2);
%             hold on;
%             plot(linspace(0,10,length(abd_f)), x, 'Color', 'black', 'LineWidth', 1);
%             if i == 8
%                 xticks(0:10)
%                 xlabel("Time(s)")
%             else
%                 xticklabels([])
%             end
%             if ii == 1
%                 ylabel(sprintf("SGC$_{%d}$", (i-1)*2+j),Interpreter="latex")
%             else
%                 ylabel(sprintf("SpSSC$_{%d}$", (i-1)*2+j),Interpreter="latex")
%             end
%             tmp = max(max(x), abs(min(x)));
%             ylim([-ceil(tmp), ceil(tmp)])
%             ax = gca;
%             p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
%             p(1) = -0.45;                % 统一的 x 偏移（越小越靠左）
%             p(2) = 0;                  % 垂直居中（0~1）
%             ax.YLabel.Position = p;
%             set(gca,'fontsize',12,'LineWidth',1.0,'TickLength',[0.01, 0.02],'Box','on')
%         end
%         if i*2 > size(sgmd_components{ii}, 2)
%             break
%         end
%     end
%     % if mod(ii, 2) == 0 || ii == 1
%     % print(sprintf('./pic/fecg_recover_method_%d.png', floor(ii/2)+1), ...
%     %     '-dpng', '-r600'); % 导出为 PNG 格式
%     % else
%     % print(sprintf('./pic/fecg_recover_method_%d_periodic.png', floor(ii/2)+1), ...
%     %     '-dpng', '-r600'); % 导出为 PNG 格式
%     % end
% end

% figure;
% delay = 8;
% sig = normalize(D.Direct_1, 'range');
% scatter(sig(1:end-delay), sig(1+delay:end))
            
% for ii = 1:7
%     fig = figure('Units','pixels', 'Position', [100, 100, 1000, 850]);
%     tlo = tiledlayout(fig,4,4, 'TileSpacing','compact', 'Padding','compact');  % 关键
%     for i = 1:4
%         for j = 1:4
%             if (i-1)*4+j > size(sgmd_components{ii}, 2)
%                 break
%             end
%             x = sgmd_components{ii}(:, (i-1)*4+j);
%             nexttile(tlo, (i-1)*4+j);
%             % xline(onset_qrs_10s, 'r--', 'LineWidth', 1.2);
%             hold on;
%             delay = 8;
%             sig = normalize(x, 'range');
%             scatter(sig(1:end-delay), sig(1+delay:end))
%             % plot(linspace(0,10,length(abd_f)), x, 'Color', 'black', 'LineWidth', 1);
%             if ii == 1
%                 ylabel(sprintf("SGC$_{%d}$", (i-1)*2+j),Interpreter="latex")
%             else
%                 ylabel(sprintf("SpSSC$_{%d}$", (i-1)*2+j),Interpreter="latex")
%             end
%             tmp = max(max(sig), abs(min(sig)));
%             ylim([0, ceil(tmp)])
%             set(gca,'fontsize',12,'LineWidth',1.0,'TickLength',[0.01, 0.02],'Box','on')
%         end
%         if i*4 > size(sgmd_components{ii}, 2)
%             break
%         end
%     end
%     if mod(ii, 2) == 0 || ii == 1
%     print(sprintf('./pic/fecg_phase_recover_method_%d.png', floor(ii/2)+1), ...
%         '-dpng', '-r600'); % 导出为 PNG 格式
%     else
%     print(sprintf('./pic/fecg_phase_recover_method_%d_periodic.png', floor(ii/2)+1), ...
%         '-dpng', '-r600'); % 导出为 PNG 格式
%     end
% end

for ii = 1:7
    RR = diff(qrs_sample);     % 样本
    T_ref = mean(RR) / fs;     % 秒
    T_k = zeros(size(sgmd_components{ii}, 2), 1);
    tmpp = T_k;
    for k = 1:size(sgmd_components{ii}, 2)-1
        xk = sgmd_components{ii}(:,k);
        xk = xk - mean(xk);
        [r,lags] = xcorr(xk,'coeff');
        lags = lags(lags>0);
        r = r(lags>0);
        [~,locs] = findpeaks(r,'MinPeakDistance',round(0.25*fs));
        T_k(k) = mean(diff(locs))/fs;
        % [~, locs] = pan_tompkin(up, fs, false);
        [~,locs] = findpeaks(xk,'MinPeakDistance',round(0.25*fs));
        tmp = eval_qrs(qrs_sample,locs,fs,10);
        tmpp(k) = tmp.F1;
    end
    % T_ref : 心跳参考周期（秒）
    % T_k   : 各分量周期（1×K 或 K×1）
    dist = abs(T_k - T_ref);     % 与心跳周期的差值
    [~, idx_sort] = sort(dist, 'ascend');
    top5_idx  = idx_sort(1:min(5,numel(dist)));        % 前 5 个分量索引
    disp([tmpp(top5_idx); mean(tmpp(top5_idx))])
    % fig = figure('Units','pixels', 'Position', [100, 100, 1200, 800]);
    % tlo = tiledlayout(fig,5,1, 'TileSpacing','compact', 'Padding','compact');  % 关键
    % for i = 1:numel(top5_idx)
    %     ax=nexttile(tlo, i);
    %     x = sgmd_components{ii}(:,top5_idx(i));
    %     % plot_envelope_spectrum(sgmd_components{ii}(:,top5_idx(i)), 1000);
    %     plot(linspace(0,10,length(abd_f)), x,'LineWidth', 0.5, 'Color','blue');
    %     hold on;
    %     xline((qrs_sample-1)/Fs, '--', 'LineWidth', .5,'Color','red');
    %     % ylim([-10,10])
    %     if i == numel(top5_idx)
    %         xlabel("time (s)")
    %     end
    %     if ii == 1
    %         ylabel(sprintf("SGC$_{%d}$", top5_idx(i)), 'Interpreter','latex','FontSize',18)
    %     else
    %         ylabel(sprintf("SGC$_{%d}$", top5_idx(i)),Interpreter="latex",FontSize=18)
    %     end
    %     if i ~= 5
    %         xticks([])
    %     end
    %     xlim([-0.15,10.15])
    %     tmp = max(max(x), abs(min(x)));
    %     ylim([-ceil(tmp), ceil(tmp)])
    %     set(ax,'fontsize',16,'LineWidth',1.0,'TickLength',[0.01, 0.02],'Box','on')
    % end
    % set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
    % if mod(ii, 2) == 0 || ii == 1
    % print(sprintf('./pic/fecg_top5_recover_method_%d.png', floor(ii/2)+1), ...
    %     '-dpng', '-r400'); % 导出为 PNG 格式
    % else
    % print(sprintf('./pic/fecg_top5_recover_method_%d_periodic.png', floor(ii/2)+1), ...
    %     '-dpng', '-r400'); % 导出为 PNG 格式
    % end
end

% delay = 8;
% color_map = lines;
% fig = figure('Units','pixels', 'Position', [100, 100, 750, 1200]);
% tlo = tiledlayout(fig,7,5, 'TileSpacing','compact', 'Padding','compact');  % 关键 
% for ii = 1:7
%     RR = diff(qrs_sample);     % 样本
%     T_ref = mean(RR) / fs;     % 秒
%     T_k = zeros(size(sgmd_components{ii}, 2), 1);
%     for k = 1:size(sgmd_components{ii}, 2)-1
%         xk = sgmd_components{ii}(:,k);
%         xk = xk - mean(xk);
%         [r,lags] = xcorr(xk,'coeff');
%         lags = lags(lags>0);
%         r = r(lags>0);
%         [~,locs] = findpeaks(r,'MinPeakDistance',round(0.25*fs));
%         T_k(k) = mean(diff(locs))/fs;
%         % [~, locs] = pan_tompkin(xk, fs, false);
%         [~,locs] = findpeaks(xk,'MinPeakDistance',round(0.25*fs));
%         tmp = eval_qrs(qrs_sample,locs,fs,10);
%         tmpp(k) = tmp.F1;
%     end
%     % T_ref : 心跳参考周期（秒）
%     % T_k   : 各分量周期（1×K 或 K×1）
%     dist = abs(T_k - T_ref);     % 与心跳周期的差值
%     [~, idx_sort] = sort(dist, 'ascend');
%     top5_idx  = idx_sort(1:min(5,numel(dist)));        % 前 5 个分量索引
%     disp(top5_idx)
%     for i = 1:numel(top5_idx)
%         ax=nexttile(tlo, (ii-1)*5+i);
%         x = sgmd_components{ii}(:,top5_idx(i));
%         sig = normalize(x, 'range');
%         scatter(sig(1:end-delay), ...
%             sig(1+delay:end), 3, color_map(k,:))
%         % plot(linspace(0,10,length(abd_f)), x,'LineWidth', 0.5, 'Color','black');
%         hold on;
%         xlabel(sprintf("SGC$_{%d}$", top5_idx(i)), 'Interpreter','latex','FontSize',18)
%         tmp = {'SGMD', 'SpSSA$_{\tau=1}$', 'CsSpSSA$_{\tau=1}$', ...
%             'SpSSA$_{\tau=2}$', 'CsSpSSA$_{\tau=2}$', ...
%             'SpSSA$_{\tau=3}$', 'CsSpSSA$_{\tau=3}$'};
%         if i == 1
%             ylabel(tmp{ii},'Interpreter','latex');
%         end
%         p = ax.XLabel.Position;   % [x y z]，归一化到当前坐标轴
%         p(1) = 0.50;                % 统一的 x 偏移（越小越靠左）
%         p(2) = -0.05;                  % 垂直居中（0~1）
%         ax.XLabel.Position = p;
%         hold on;
%         xlim([0,1])
%         ylim([0,1])
%         xticks([0,1])
%         yticks([0,1])
%         set(ax,'fontsize',10,'LineWidth',1.0,'TickLength',[0.01, 0.02],'Box','on')
%     end
% end
% set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
% print('./pic/fecg_top5_phase_recover_methods.png', '-dpng', '-r400'); % 导出为 PNG 格式

end

%% 1) 读 CSV
function D = read_fecg_csv(csvFile)
% 该CSV第1行是单位行（uV），第2行起是数据
opts = detectImportOptions(csvFile, 'Delimiter', ',');
T = readtable(csvFile, opts);

% 第一行是单位，剔除
% T = T(2:end, :);

% 列名可能带引号/空格，做一次鲁棒映射
vars = T.Properties.VariableNames;

% 辅助函数：按包含字段找列
col = @(name) find(contains(vars, name, 'IgnoreCase', true), 1, 'first');

iTime   = col('Time');
iDirect = col('Direct');
iA1     = col('Abdomen_1');
iA2     = col('Abdomen_2');
iA3     = col('Abdomen_3');
iA4     = col('Abdomen_4');

% 转数值（table里是字符串时需要 str2double）
D.Time = T{:, iTime}; % 时间戳保留（通常不用）
D.Direct_1  = toNum(T{:, iDirect});
D.Abdomen_1 = toNum(T{:, iA1});
D.Abdomen_2 = toNum(T{:, iA2});
D.Abdomen_3 = toNum(T{:, iA3});
D.Abdomen_4 = toNum(T{:, iA4});
end

function x = toNum(col)
if iscell(col) || isstring(col)
    x = str2double(string(col));
else
    x = double(col);
end
x = x(:);
end

%% 2) 预处理
function y = preprocess_ecg(x, Fs, useNotch, notchHz)
x = x(:);

% 去均值 + 去趋势
x = x - mean(x, 'omitnan');
x = detrend(x);

% 带通 1-100 Hz（可按你论文设定调整）
bp = designfilt('bandpassiir', ...
    'FilterOrder', 6, ...
    'HalfPowerFrequency1', 1, ...
    'HalfPowerFrequency2', 100, ...
    'SampleRate', Fs);
y = filtfilt(bp, x);

% 陷波（50 or 60 Hz）
if useNotch
    wo = notchHz/(Fs/2);
    bw = wo/35; % 经验带宽，想更窄可调大分母
    [b,a] = iirnotch(wo, bw);
    y = filtfilt(b, a, y);
end

% 归一化（可选：让相关/阈值更稳定）
y = y ./ (std(y) + eps);
end

function [xSUM, w, info] = build_weighted_sum_GSVD(D, fs)
%BUILD_WEIGHTED_SUM_GSVD
%   Construct weighted abdominal ECG sum to enhance MECG
%   using GSVD-based weighting (Jezewski-style).
%
% INPUT:
%   D  - struct with fields:
%        D.Abdomen_1, D.Abdomen_2, D.Abdomen_3, D.Abdomen_4
%   fs - sampling frequency (Hz), e.g. 1000
%
% OUTPUT:
%   xSUM - weighted composite signal (MECG enhanced)
%   w    - weight vector (4x1)
%   info - structure with GSVD diagnostics

% -------- 1. Stack abdominal signals --------
X = [ ...
    D.Abdomen_1(:), ...
    D.Abdomen_2(:), ...
    D.Abdomen_3(:), ...
    D.Abdomen_4(:) ...
    ];

% Remove mean to avoid DC dominance
X = X - mean(X, 1);

% -------- 2. Extract maternal-dominant component --------
% Maternal ECG approx. 5–40 Hz
[bM, aM] = butter(4, [5 50] / (fs/2), 'bandpass');
XM = filtfilt(bM, aM, X);

% -------- 3. Residual = non-maternal subspace --------
XN = X - XM;

% -------- 4. GSVD on covariance matrices (FIXED) --------
% Covariance matrices (4 × 4)
RM = cov(XM);   % maternal-dominant covariance
RN = cov(XN);   % non-maternal covariance

% GSVD on covariance level
[~, V, ~, C, S] = gsvd(RM, RN);

c = diag(C);
s = diag(S);
ratio = c ./ s;

[~, idx] = max(ratio);
w = V(:, idx);

% Ensure column vector and normalize
w = w(:);
w = w / norm(w);

% -------- 5. Weighted sum signal --------
xSUM = X * w;

% -------- 6. Output diagnostics --------
info = struct();
info.weight_ratio = ratio;
info.selected_index = idx;
info.weights = w;
info.maternal_energy = var(XM * w);
info.residual_energy = var(XN * w);

end


function S = eval_qrs(ref, det, fs, tol_ms)

ref = sort(ref(:));
det = sort(det(:));
tol = round(fs * tol_ms / 1000);

i = 1; j = 1;
TP = 0; FP = 0; FN = 0;
err = [];   % matched timing errors (samples)

while i <= numel(ref) && j <= numel(det)
    d = det(j) - ref(i);

    if abs(d) <= tol
        TP = TP + 1;
        err(end+1,1) = d; %#ok<AGROW>
        i = i + 1;
        j = j + 1;
    elseif det(j) < ref(i) - tol
        FP = FP + 1;
        j = j + 1;
    else
        FN = FN + 1;
        i = i + 1;
    end
end

FN = FN + (numel(ref) - i + 1);
FP = FP + (numel(det) - j + 1);

Se  = TP / max(TP+FN, 1);
PPV = TP / max(TP+FP, 1);
F1  = 2*TP / max(2*TP + FP + FN, 1);

err_ms = 1000*err/fs;
S.TP = TP; S.FP = FP; S.FN = FN;
S.Se = Se; S.PPV = PPV; S.F1 = F1;

if ~isempty(err_ms)
    S.MedAE_ms = median(abs(err_ms));
    S.MAE_ms   = mean(abs(err_ms));
    S.Bias_ms  = mean(err_ms);
    S.P95_ms   = prctile(abs(err_ms), 95);
else
    S.MedAE_ms = NaN; S.MAE_ms = NaN; S.Bias_ms = NaN; S.P95_ms = NaN;
end
end
