
clc
clear
close all

% signal 1
fs = 5120;
time = 1;
t = 0:1/fs:time-1/fs;
x = zeros(fs, 3);
x(:, 1) = sin(120*pi*t);
x(:, 2) = sin(60*pi*t).*sin(2*pi*t);
x(:, 3) = 0.5*cos(10*pi*t).*(1+0.5*sin(50*pi*t));
% x(:, 4) = 0.5*t+0.1*t.*t;
sig = sum(x, 2);

%%
tau = 1;
n_cluster = 3;
threshold_corr = 0.85;  % 0.85
threshold_nmse = 0.01;  % 0.01
threshold_svd = 1;  % 1
SGCs = cell(2, 1);
SGCs{1}=SGMD(sig',fs,1,threshold_corr,threshold_nmse);
% SGCs{1}=DRSGMD(sig,t,fs,1,2,threshold_corr,threshold_nmse);
periodic_extension = false;
SGCs{2} = SpSSA(sig, fs, 1, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
opts.threshold_nmse=0.05; SGCs{2} = SSD(sig, fs, 2, opts);
% SGCs{4} = SpSSA(sig, fs, 2, n_cluster, ...
%     threshold_corr, threshold_nmse, threshold_svd, ...
%     'kmeans', periodic_extension);
% SGCs{6} = SpSSA(sig, fs, 3, n_cluster, ...
%     threshold_corr, threshold_nmse, threshold_svd, ...
%     'kmeans', periodic_extension);
% periodic_extension = true;
% SGCs{3} = SpSSA(sig, fs, 1, n_cluster, ...
%     threshold_corr, threshold_nmse, threshold_svd, ...
%     'kmeans', periodic_extension);
% SGCs{5} = SpSSA(sig, fs, 2, n_cluster, ...
%     threshold_corr, threshold_nmse, threshold_svd, ...
%     'kmeans', periodic_extension);
% SGCs{7} = SpSSA(sig, fs, 3, n_cluster, ...
%     threshold_corr, threshold_nmse, threshold_svd, ...
%     'kmeans', periodic_extension)
for i = 1:numel(SGCs)
    if size(x, 2) ~= size(SGCs{i}, 2)-1
        continue
    end
    idx = [alignSignalsByCorr(x, SGCs{i}(:, 1:end-1)), size(SGCs{i}, 2)];
    SGCs{i} = SGCs{i}(:, idx);
end

relative_fun = @(x, y) norm(y-x, 2) / norm(x, 2);
rmse_fun = @(x, y) sqrt(mean((y - x).^2));
% 全量（把输入展平成向量后计算一个标量 MAE）
nmae_fun = @(x, y) mean(abs(y(:) - x(:))) / max(eps, (max(x(:)) - min(x(:))));
rMat = zeros(size(x, 2), numel(SGCs));
nMat = zeros(size(x, 2), numel(SGCs));
for i = 1:size(x, 2)
    for j = 1:numel(SGCs)
        rMat(i, j) = relative_fun(x(:, i), SGCs{j}(:, i));
        fprintf("%.4f\t", rMat(i, j));
    end
    disp(" ")
end
fprintf("%.4f\t", mean(rMat, 1));
disp(" ")
disp(" ")
for i = 1:size(x, 2)
    for j = 1:numel(SGCs)
        nMat(i, j) = nmae_fun(x(:, i), SGCs{j}(:, i));
        fprintf("%.4f\t", nMat(i, j));
    end
    disp(" ")
end
fprintf("%.4f\t", mean(nMat, 1));
disp(" ")
disp(" ")
for i = 1:size(x, 2)
    for j = 1:numel(SGCs)
        fprintf("%.4f\t", corr(x(:, i), SGCs{j}(:, i)));
    end
    disp(" ")
end

%%
color_map = lines;  % colormap used for lines and markers
labels = {'$x_1(t)$', '$x_2(t)$', '$x_3(t)$', '$x_4(t)$', '$n(t)$', '$s(t)$'};
fig = figure('Units','pixels', 'Position', [100, 100, 1200, 800]);
tlo = tiledlayout(fig,5,1, 'TileSpacing','compact', 'Padding','compact');
for i = 1:size(x, 2)
    ax = nexttile(tlo,i);
    % plot_envelope_spectrum(x(:,i), fs);
    plot(t, x(:, i),'LineWidth',1.0,'Color','blue');
    xlim([0,time]);
    ylabel(labels{i}, 'Interpreter', 'latex','FontSize',18)
    set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
    p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
    p(1) = -0.04;                % 统一的 x 偏移（越小越靠左）
    p(2) = 0;                  % 垂直居中（0~1）
    ax.YLabel.Position = p;
end
ax = nexttile(tlo,5);
plot(t, sig,'LineWidth',1.0,'Color','blue');
xlim([0,time]);
ylabel(labels{6}, 'Interpreter', 'latex','FontSize',18)
set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
p(1) = -0.04;                % 统一的 x 偏移（越小越靠左）
p(2) = 0.5;                  % 垂直居中（0~1）
ax.YLabel.Position = p;
xlabel('time(s)', 'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
% print('./pic/synthetic_signal.png', '-dpng', '-r400'); % 导出为 PNG 格式

for j = 1:numel(SGCs)
    b = size(SGCs{j}, 2);
    fig=figure('Units','pixels', 'Position', [100, 100, 1200, 800]);
    tlo = tiledlayout(fig,b,1, 'TileSpacing','compact', 'Padding','compact');  % 关键
    for i = 1:b
        ax=nexttile(tlo, i);
        % plot_envelope_spectrum(SGCs{j}(:,i),fs);
        plot(t, SGCs{j}(:,i),'LineWidth',1.0,'Color','blue');
        xlim([0,time]);
        ylabel(labels{i}, 'Interpreter', 'latex','FontSize',18)
        set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
        p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
        p(1) = -0.04;                % 统一的 x 偏移（越小越靠左）
        p(2) = sum(ax.YLim) / 2;                  % 垂直居中（0~1）
        ax.YLabel.Position = p;
    end
    xlabel('time(s)', 'FontSize',16)
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
    % if mod(j, 2) == 0 || j == 1
    % print(sprintf('./pic/recover_method_%d.png', floor(j/2)+1), ...
    %     '-dpng', '-r400'); % 导出为 PNG 格式
    % else
    % print(sprintf('./pic/recover_method_%d_periodic.png', floor(j/2)+1), ...
    %     '-dpng', '-r400'); % 导出为 PNG 格式
    % end
end
