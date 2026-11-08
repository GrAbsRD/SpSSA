
clc
clear
close all

% signal 1
fs = 5120;
time = 1;
t = 0:1/fs:time-1/fs;
x = zeros(fs, 3);
x(:, 1) = 2*(1+0.5*(sin(2*pi*t))).*sin(60*pi*t);
x(:, 2) = sin(120*pi*t);
x(:, 3) = 0.5*cos(10*pi*t);
sig = sum(x, 2);

%% signal 2
% fs = 1000;          % sampling rate [Hz]
% time = 1;
% t = 0:1/fs:time-1/fs;
% x = zeros(fs, 3);
% % Components
% x(:, 1) = sin(2*pi*5*t);
% % small "carrier" ripple at 15 Hz with amplitude 0.1
% x(:, 2) = 0.1*sin(2*pi*15*t);
% % high-frequency intermittent (75 Hz) added only on the second half.
% gate = (t >= 0.5);                          % only second half of the record
% x(:, 3) = gate .* sin(2*pi*75*t);    % 75 Hz, unit-amp bursts
% sig = sum(x, 2);

%%
tau = 1;
n_cluster = 3;
threshold_corr = 0.85;  % 0.85
threshold_nmse = 0.01;  % 0.01
threshold_svd = 1;  % 1
SGCs = cell(7, 1);
SGCs{1}=SGMD(sig',fs,1,threshold_corr,threshold_nmse);
% SGCs{1}=vmd(sig','NumIMFs',n_cluster);
periodic_extension = false;
% SGCs{2} = SGMD_SVD(sig, fs, 1, n_cluster, ...
%     threshold_corr, threshold_nmse, ...
%     'kmeans', periodic_extension);
SGCs{2} = SpSSA(sig, fs, 4, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
% opts.threshold_nmse=0.001; SGCs{2} = SSD(sig, fs, 2, opts);
% tmp = JMD(sig,5000,50,0.03,0.45,0.8,4,2,1e-6);
SGCs{4} = SpSSA(sig, fs, 5, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
SGCs{6} = SpSSA(sig, fs, 10, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
periodic_extension = true;
SGCs{3} = SpSSA(sig, fs, 4, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
SGCs{5} = SpSSA(sig, fs, 5, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
SGCs{7} = SpSSA(sig, fs, 10, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension)
for i = 1:numel(SGCs)
    if size(x, 2) ~= size(SGCs{i}, 2)-1
        continue
    end
    idx = [alignSignalsByCorr(x, SGCs{i}(:, 1:end-1)), size(SGCs{i}, 2)];
    SGCs{i} = SGCs{i}(:, idx);
end

relative_fun = @(x, y) norm(y-x, 2) / norm(x, 2);
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
% component_idx = 3;
% labels = {'SpSSA$_{\tau=1}$', 'SpSSA$_{\tau=2}$', 'SpSSA$_{\tau=3}$'};
% 
% fig=figure('Units','pixels', 'Position', [100, 100, 1200, 870]);
% tlo = tiledlayout(fig,3,1, 'TileSpacing','compact', 'Padding','compact');  % 关键
% for i = 1:3
%     ax=nexttile(tlo, i);
%     plot(t, x(:, component_idx),'--','LineWidth',1.0,'Color','black');
%     hold on;
%     plot(t, SGCs{i*2}(:,component_idx),'LineWidth',1.0,'Color','blue');
%     hold on;
%     plot(t, SGCs{i*2+1}(:,component_idx),'LineWidth',1.0,'Color','red');
%     hold on;
%     plot(t, x(:, component_idx),'--','LineWidth',1.0,'Color','black');
%     xlim([0,time]);
%     ylim([-0.7,0.7])
%     yticks([-0.5,0,0.5])
%     yticklabels({'-0.5', '0', '0.5'})
%     ylabel(labels{i}, 'Interpreter', 'latex')
%     set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
%     p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
%     p(1) = -0.04;                % 统一的 x 偏移（越小越靠左）
%     p(2) = sum(ax.YLim) / 2;                  % 垂直居中（0~1）
%     ax.YLabel.Position = p;
%     if i == 3
%         legend({sprintf('$x_%d(t)$', component_idx), ...
%             sprintf('SpSSC$_{%d}$', component_idx), ...
%             sprintf('SpSSC$_{%d{\\rm -ext}}$', component_idx)}, ...
%             'Interpreter','latex','NumColumns',3, ...
%             'Location','southoutside','Box','off','FontSize',16)
%     end
% end
% xlabel('time(s)', 'FontSize',16)
% set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
% print(sprintf('./pic/comparision_x%d.png', component_idx), '-dpng', '-r400'); % 导出为 PNG 格式



%%
color_map = lines;  % colormap used for lines and markers
labels = {'$x_1(t)$', '$x_2(t)$', '$x_3(t)$', 'noise', '$x(t)$'};
fig = figure('Units','pixels', 'Position', [100, 100, 1200, 800]);
tlo = tiledlayout(fig,4,1, 'TileSpacing','compact', 'Padding','compact');
for i = 1:size(x, 2)
    ax = nexttile(tlo,i);
    plot(t, x(:, i),'LineWidth',1.0,'Color','blue');
    xlim([0,time]);
    ylabel(labels{i}, 'Interpreter', 'latex','FontSize',18)
    set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
    p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
    p(1) = -0.04;                % 统一的 x 偏移（越小越靠左）
    p(2) = 0;                  % 垂直居中（0~1）
    ax.YLabel.Position = p;
end
ax = nexttile(tlo,4);
plot(t, sig,'LineWidth',1.0,'Color','blue');
xlim([0,time]);
ylabel(labels{5}, 'Interpreter', 'latex','FontSize',18)
set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
p(1) = -0.04;                % 统一的 x 偏移（越小越靠左）
p(2) = 0.5;                  % 垂直居中（0~1）
ax.YLabel.Position = p;
xlabel('time(s)', 'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
print('./pic/synthetic_signal.png', '-dpng', '-r400'); % 导出为 PNG 格式

for j = [1 2 4 6]
    b = size(SGCs{j}, 2);
    fig=figure('Units','pixels', 'Position', [100, 100, 1200, 850]);
    tlo = tiledlayout(fig,b,1, 'TileSpacing','compact', 'Padding','compact');  % 关键
    for i = 1:b
        ax=nexttile(tlo, i);
        if j == 1
            plot(t, SGCs{j}(:,i),'LineWidth',1.0,'Color','red');
            hold on;
        else
            plot(t, SGCs{j}(:,i),'LineWidth',1.0,'Color','blue');
            hold on;
            plot(t, SGCs{j+1}(:,i),'LineWidth',1.0,'Color','red');
            hold on;
        end
        if i ~= b
            plot(t, x(:, i),'--','LineWidth',1.0,'Color','black');
            hold on;
        else
            plot(t, zeros(fs,1),'--','LineWidth',1.0,'Color','black');
            hold on;
        end
        xlim([0,time]);
        if i == 3
            ylim([min(ax.YLim(1), -0.5), max(ax.YLim(2), 0.5)]);
            yticks([-0.5,0,0.5])
            yticklabels({'-0.5', '0', '0.5'})
        end
        if i == 4 && j == 1
            ax.YAxis.ExponentMode = 'manual';
            ax.YAxis.Exponent = -2;   % 相当于整体乘 10^3
            legend({'SGCs','Original signal'}, ...
                'Interpreter','latex','NumColumns',2, ...
                'Location','southoutside','Box','off','FontSize',16)
        end
        ylabel(labels{i}, 'Interpreter', 'latex','FontSize',18)
        set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
        p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
        p(1) = -0.04;                % 统一的 x 偏移（越小越靠左）
        p(2) = sum(ax.YLim) / 2;                  % 垂直居中（0~1）
        ax.YLabel.Position = p;
        if i == 4 && j ~= 1
            legend({'SGCs','CsSGCs','Original signal'}, ...
                'Interpreter','latex','NumColumns',3, ...
                'Location','southoutside','Box','off','FontSize',16)
        end
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
    print(sprintf('./pic/recover_method_%d_all.png', floor(j/2)+1), ...
        '-dpng', '-r400'); % 导出为 PNG 格式
end
