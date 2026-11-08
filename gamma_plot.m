%% Grouped horizontal bar chart with value labels + bracketed mean
clear; clc; close all;

%% ====== Data from the LaTeX table ======
components = {'x_1(t)','x_2(t)'};

results = [0.2450    0.4827
    0.2450    0.4827
    0.2450    0.4827
    0.2450    0.4827
    0.2450    0.0660
    0.2450    0.0660
    0.2450    0.0660
    0.2450    0.0660
    0.2450    0.0660
    0.2450    0.0660
    0.2450    0.0660
    0.2450    0.0660
    0.2450    0.0660
    0.2450    0.1936
    0.2450    0.1936
    0.2450    0.3675
    0.2450    0.3675
    0.2450    0.4659
    0.2450    0.4969
    0.2450    0.5337
    0.2450    0.5662];

% Averages (use table values or compute)
avgVals = [0.3639
    0.3639
    0.3639
    0.3639
    0.1555
    0.1555
    0.1555
    0.1555
    0.1555
    0.1555
    0.1555
    0.1555
    0.1555
    0.2193
    0.2193
    0.3062
    0.3062
    0.3554
    0.3709
    0.3894
    0.4056];
filename = 'relative_error'; xLabel = 'Relative error'; plength=900; pwidth=1000;

%%
fig = figure('Units','pixels', 'Position', [100, 100, 850, 500]);

gamma = 0.6:0.01:0.8;

color_map = lines;
% color_map = {'black','blue','blue','blue','red','red','red'}
markers = {'o','s','^','d','v','>','p'};
for i = 1:2
    idx = 1:21;
    if i == 2
        lw = 1.3;
        ms = 7;
    else
        lw = 1;
        ms = 6;
    end
    plot(gamma(idx), results(idx,i), ...
        '-','LineWidth',lw, ...
        'Color',color_map(i,:), ...
        'Marker',markers{i}, ...
        'MarkerFaceColor',[1 1 1], ...   % 白心
        'MarkerEdgeColor',color_map(i,:), ...
        'MarkerSize',ms);
    hold on
end

i = 3;
plot(gamma(idx), avgVals(idx), ...
        '--','LineWidth',lw, ...
        'Color',color_map(i,:), ...
        'Marker',markers{i}, ...
        'MarkerFaceColor',[1 1 1], ...   % 白心
        'MarkerEdgeColor',color_map(i,:), ...
        'MarkerSize',6);

% 每一条曲线对应一个名字
lgd_names = { ...
    '$x_{1}(t)$', ...            % i = 1
    '$x_{2}(t)$', ... % i = 2
    'avg' ...
};
xlim([0.59, 0.81])
xticks(gamma)
% labels = string(gamma);
% labels(mod(dd,50)~=0) = "";
% xticklabels(labels)
xtickangle(90)
xlabel("$\gamma$", 'Interpreter','latex')
ylabel("Relative error")
legend(lgd_names, 'Interpreter','latex','Location', 'southeast','FontSize',16);
set(gca,'fontsize',16,'LineWidth',1.0,'TickLength',[0.01, 0.02],'Box','on')
set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
print('./pic/gamma_relative_error.png', '-dpng', '-r400');