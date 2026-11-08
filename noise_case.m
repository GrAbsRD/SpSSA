
clc
clear
close all

fs = 5120;
time = 1;
t = 0:1/fs:time-1/fs;
x = zeros(fs, 2);
x(:, 1) = 2*(1+0.5*(sin(2*pi*t))).*sin(40*pi*t);
x(:, 2) = 1.2*cos(25*t);

rng(56); % 固定种子
% 总信号（不加噪声）
sig_clean = sum(x, 2);
% 目标 SNR（dB）
target_snr_db = 5;
% 信号标准差（总体标准差，ddof=0）
sigma_x = std(sig_clean);
% 计算所需噪声标准差（见上式）
sigma_n = sigma_x / sqrt(10^(target_snr_db / 10));
% 生成噪声并加到信号
noise = randn(length(t), 1) * sigma_n;
sig_noisy = sig_clean + noise;


%%
tau = 1;
n_cluster = 2;
threshold_corr = 0.7;
threshold_nmse = 0.3;
threshold_svd = 0.75;

SGCs = cell(7, 1);
SGCs{1}=SGMD(sig_noisy',fs,1,threshold_corr,threshold_nmse);
periodic_extension = false;
SGCs{2} = SpSSA(sig_noisy, fs, 4, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
% opts.threshold_nmse=threshold_nmse; SGCs{2} = SSD(sig_noisy, fs, 3, opts)
SGCs{4} = SpSSA(sig_noisy, fs, 5, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
SGCs{6} = SpSSA(sig_noisy, fs, 10, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
periodic_extension = true;
SGCs{3} = SpSSA(sig_noisy, fs, 4, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
SGCs{5} = SpSSA(sig_noisy, fs, 5, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);
SGCs{7} = SpSSA(sig_noisy, fs, 10, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);

% SGCs{3} = movmean(SGCs{3}, 10);
% SGCs{4} = movmean(SGCs{4}, 10);

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
% component_idx = 1;
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
%     if component_idx == 2
%     ylim([-1.5,1.5])
%     yticks([-1.0,0,1.0])
%     yticklabels({'-1.0', '0', '1.0'})
%     end
%     ylabel(labels{i}, 'Interpreter', 'latex')
%     set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
%     p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
%     p(1) = -0.04;                % 统一的 x 偏移（越小越靠左）
%     p(2) = sum(ax.YLim) / 2;                  % 垂直居中（0~1）
%     ax.YLabel.Position = p;
%     if i == 3
%         legend({sprintf('$x_%d(t)$', component_idx), ...
%             sprintf('SpSSC$_{%d}$', component_idx), ...
%             sprintf('SpSSC$_{%d{\\rm -ext}}$', component_idx)},...
%         'Interpreter','latex','NumColumns',3, ...
%         'Location','southoutside','Box','off','FontSize',16)
%     end
% end
% xlabel('time(s)', 'FontSize',16)
% set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
% print(sprintf('./pic/comparision_x%d_noise.png', component_idx), '-dpng', '-r400'); % 导出为 PNG 格式


%%
color_map = lines;  % colormap used for lines and markers
labels = {'$x_1(t)$', '$x_2(t)$', '$noise$', '$x_o(t)$'};
fig = figure('Units','pixels', 'Position', [100, 100, 1200, 800]);
tlo = tiledlayout(fig,3,1, 'TileSpacing','compact', 'Padding','compact');
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
ax = nexttile(tlo,3);
plot(t, sig_noisy,'LineWidth',1.0,'Color','blue');
xlim([0,time]);
ylabel(labels{4}, 'Interpreter', 'latex','FontSize',18)
set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
p(1) = -0.04;                % 统一的 x 偏移（越小越靠左）
p(2) = 0.5;                  % 垂直居中（0~1）
ax.YLabel.Position = p;
xlabel('time(s)', 'FontSize',16)
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
print('./pic/synthetic_signal_noise.png', '-dpng', '-r400'); % 导出为 PNG 格式

for j = [1 2 4 6]
    b = size(SGCs{j}, 2) - 1;
    fig=figure('Units','pixels', 'Position', [100, 100, 1200, 550]);
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
        plot(t, x(:, i),'--','LineWidth',1.0,'Color','black');
        hold on;
        xlim([0,time]);
        if i == 2
            ylim([min(ax.YLim(1), -1.5), max(ax.YLim(2), 1.5)]);
            yticks([-1.5,0,1.5])
            yticklabels({'-1.5', '0', '1.5'})
        end
        ylabel(labels{i}, 'Interpreter', 'latex','FontSize',18)
        set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
        p = ax.YLabel.Position;   % [x y z]，归一化到当前坐标轴
        p(1) = -0.04;                % 统一的 x 偏移（越小越靠左）
        p(2) = sum(ax.YLim) / 2;                  % 垂直居中（0~1）
        ax.YLabel.Position = p;
        if i == 2 && j == 1
            legend({'SGCs','Original signal'}, ...
                'Interpreter','latex','NumColumns',2, ...
                'Location','southoutside','Box','off','FontSize',16)
        end
        if i == 2 && j ~= 1
            legend({'SGCs','CsSGCs','Original signal'}, ...
                'Interpreter','latex','NumColumns',3, ...
                'Location','southoutside','Box','off','FontSize',16)
        end
    end
    xlabel('time(s)', 'FontSize',16)
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
    % if mod(j, 2) == 0 || j == 1
    % print(sprintf('./pic/noise_recover_method_%d.png', floor(j/2)+1), ...
    %     '-dpng', '-r400'); % 导出为 PNG 格式
    % else
    % print(sprintf('./pic/noise_recover_method_%d_periodic.png', floor(j/2)+1), ...
    %     '-dpng', '-r400'); % 导出为 PNG 格式
    % end
    print(sprintf('./pic/noise_recover_method_%d_all.png', floor(j/2)+1), ...
        '-dpng', '-r400'); % 导出为 PNG 格式
end

%%
for j = [2, 4]
    b = size(SGCs{j}, 2) - 1;
    fig=figure('Units','pixels', 'Position', [100, 100, 1200, 800]);
    tlo = tiledlayout(fig,b,1, 'TileSpacing','compact', 'Padding','compact');  % 关键
    for i = 1:b    
        ax=nexttile(tlo, i);  
        plot(t, SGCs{j}(:,i),'LineWidth',1.0,'Color','blue');  
        xlim([0.7,0.8]);
            ylim([min(ax.YLim(1), 0), max(ax.YLim(2), 1.5)]);
            yticks([0,1.5])
            yticklabels({'0', '1.5'})
        ylabel(labels{i}, 'Interpreter', 'latex','FontSize',18)
        set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
    end
    xlabel('time(s)', 'FontSize',16)
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
    print(sprintf('./pic/noise_detail_method_%d.png', floor(j/2)+1), '-dpng', '-r400'); % 导出为 PNG 格式
end
% for j = 5:2:7
%     b = size(SGCs{j}, 2) - 1;
%     fig=figure('Units','pixels', 'Position', [100, 100, 1200, 800]);
%     tlo = tiledlayout(fig,b,1, 'TileSpacing','compact', 'Padding','compact');  % 关键
%     for i = 1:b    
%         ax=nexttile(tlo, i);  
%         plot(t, SGCs{j}(:,i),'LineWidth',1.0,'Color','blue');  
%         xlim([0.7,0.8]);
%             ylim([min(ax.YLim(1), 0), max(ax.YLim(2), 1.5)]);
%             yticks([0,1.5])
%             yticklabels({'0', '1.5'})
%         ylabel(labels{i}, 'Interpreter', 'latex','FontSize',18)
%         set(ax,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])
%     end
%     xlabel('time(s)', 'FontSize',16)
%     set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
%     print(sprintf('./pic/noise_detail_method_%d_periodic.png', floor(j/2)+1), '-dpng', '-r400'); % 导出为 PNG 格式
% end


%%
tau = 1;
n_cluster = 2;
threshold_corr = 0.7;
threshold_nmse = 0.3;

tmp = 0.6:0.01:0.8;

SGCs = cell(length(tmp), 1);
periodic_extension = true;

for i = 1:length(tmp)

threshold_svd = tmp(i);
SGCs{i}=SpSSA(sig_noisy, fs, 10, n_cluster, ...
    threshold_corr, threshold_nmse, threshold_svd, ...
    'kmeans', periodic_extension);

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

