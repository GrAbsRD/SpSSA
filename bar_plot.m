%% Grouped horizontal bar chart with value labels + bracketed mean
clear; clc; close all;

%% ====== Data from the LaTeX table ======
components = {'x_1(t)','x_2(t)','x_3(t)'};

% SGMD = [0.0382; 0.0669; 0.1226];
% Sp1_non = [0.0382; 0.0669; 0.1209];
% Sp1_per = [0.0245; 0.0503; 0.0862];
% Sp2_non = [0.0322; 0.0354; 0.1224];
% Sp2_per = [0.0149; 0.0284; 0.0330];
% Sp3_non = [0.0532; 0.0168; 0.2232];
% Sp3_per = [0.0166; 0.0132; 0.0650];
% Sp4_non = [0.0781; 0.0091; 0.3279];
% Sp4_per = [0.0270; 0.0042; 0.1130];
% Sp5_non = [0.1032; 0.0047; 0.4317];
% Sp5_per = [0.0410; 0.0009; 0.1723];
% % Averages (use table values or compute)
% avgVals = [0.0759, 0.0754, 0.0537, 0.0633, 0.0254, 0.0977, 0.0316, 0.1384, 0.0481, 0.1799, 0.0714];
% filename = 'relative_error'; xLabel = 'Relative error'; plength=900; pwidth=1000;

% NMAE (from Table tab:MAE)
SGMD    = [0.0061; 0.0174; 0.0305];
Sp1_non = [0.0061; 0.0174; 0.0303];
Sp1_per = [0.0051; 0.0156; 0.0273];
Sp2_non = [0.0047; 0.0094; 0.0210];
Sp2_per = [0.0031; 0.0089; 0.0096];
Sp3_non = [0.0074; 0.0043; 0.0435];
Sp3_per = [0.0033; 0.0041; 0.0187];
Sp4_non = [0.0118; 0.0018; 0.0703];
Sp4_per = [0.0055; 0.0159; 0.0023];
Sp5_non = [0.0165; 0.0006; 0.0978];
Sp5_per = [0.0083; 0.0002; 0.0494];
% Averages (table values)
avgVals = [0.0180, 0.0179, 0.0160, 0.0117, 0.0072, 0.0184, 0.0087, 0.0280, 0.0079, 0.0383, 0.0193];
filename = 'nMAE'; xLabel = 'nMAE'; plength=900; pwidth=1000;

% Relative reconstruction error (noise case) from Table tab:relative_error_noise
% SGMD    = [0.0886; 0.1503];
% Sp1_non = [0.0886; 0.1503];
% Sp1_per = [0.0692; 0.1132];
% Sp2_non = [0.0757; 0.0794];
% Sp2_per = [0.0408; 0.0693];
% Sp3_non = [0.1198; 0.0866];
% Sp3_per = [0.0483; 0.0727];
% % Averages (table values)
% avgVals = [0.1194, 0.1194, 0.0912, 0.0776, 0.0550, 0.1032, 0.0605];
% filename = 'relative_error_noise'; xLabel = 'Relative error'; plength=900; pwidth=800;

% NMAE (noise case) from Table tab:MAE_noise
% SGMD    = [0.0167; 0.0373];
% Sp1_non = [0.0167; 0.0373];
% Sp1_per = [0.0144; 0.0328];
% Sp2_non = [0.0126; 0.0196];
% Sp2_per = [0.0080; 0.0180];
% Sp3_non = [0.0194; 0.0227];
% Sp3_per = [0.0098; 0.0204];
% % Averages (table values)
% avgVals = [0.0270, 0.0270, 0.0236, 0.0161, 0.0130, 0.0210, 0.0151];
% filename = 'nMAE_noise'; xLabel = 'nMAE'; plength=900; pwidth=800;

%% ====== Layout settings ======
methodGroups = { ...
    'SGMD', ...
    'SpSSA_{\tau=1}', 'CsSpSSA_{\tau=1}', ...
    'SpSSA_{\tau=2}', 'CsSpSSA_{\tau=2}', ...
    'SpSSA_{\tau=3}', 'CsSpSSA_{\tau=3}', ...
    'SpSSA_{\tau=4}', 'CsSpSSA_{\tau=4}', ...
    'SpSSA_{\tau=5}', 'CsSpSSA_{\tau=5}'};

dataCell = {SGMD, Sp1_non, Sp1_per, Sp2_non, Sp2_per, Sp3_non, Sp3_per, Sp4_non, Sp4_per, Sp5_non, Sp5_per};

nGroups = numel(methodGroups);
nComp   = numel(SGMD); % 3 bars per group

gapInGroup = 0.0;   % spacing between bars in the same group
gapBetween = 0.36;   % spacing between method groups
barH       = 0.32;   % bar thickness (height in y direction)

%% ====== Build positions for a horizontal barh-like plot (manual) ======
% We'll draw each bar using rectangle() for full control of y-positions.
yPos = zeros(nGroups, nComp);
yTickPos = zeros(nGroups,1);

y0 = 0; % start
for g = 1:nGroups
    % bars in this group: stacked vertically (y changes), each bar is one component
    % top to bottom: x1, x2, x3 (you can flip if desired)
    local = y0 + (0:nComp-1)*(barH + gapInGroup);
    yPos(g,:) = local;
    yTickPos(g) = mean(local); % group center for labeling
    y0 = local(end) + barH + gapBetween;
end

%% ====== Plot ======
figure('Units','pixels', 'Position', [100, 100, plength, pwidth]);
hold on;

% Choose x-limits a bit larger than max for labels/brackets
% ===== Robust xMax =====
allVals = vertcat(dataCell{:});   % 变成 (3*nGroups) x 1
allVals = double(allVals(:));     % 强制数值列向量

xMax = max(allVals, [], 'omitnan');  % 忽略 NaN（需要较新版本）
xMax = xMax * 1.35;

% 保护：避免 NaN/Inf/非正数导致 xlim 报错
if isempty(xMax) || ~isfinite(xMax) || xMax <= 0
    xMax = 1; % fallback
end

color_map = lines;
% Draw bars
for g = 1:nGroups
    vals = dataCell{g};
    for k = 1:nComp
        % rectangle at [x y w h] with x from 0 to value
        x = 0;
        y = yPos(g,k);
        w = vals(k);
        h = barH;
        rectangle('Position',[x y w h], ...
                  'FaceColor',(1-color_map(k,:))/4.5+color_map(k,:), ...
                  'LineWidth', 1., ...
                  'HandleVisibility','off');

        % value text to the right of bar
        txt = sprintf('%.4f', w);
        
        xTxt = double(w + 0.01*xMax);
        yTxt = double(y + h/2);
        
        % 强制标量（防止 w 或 y 意外成向量）
        xTxt = xTxt(1);
        yTxt = yTxt(1);
        
        text(xTxt, yTxt, txt, ...
             'VerticalAlignment','middle', ...
             'HorizontalAlignment','left', ...
             'FontSize',10);
    end
end

% y-axis labels = method groups
set(gca,'YTick',yTickPos,'YTickLabel',methodGroups,'FontSize',16);
% ylabel('Method'); 
xlabel(xLabel);
xlim([0 xMax]);

% Put x1/x2/x3 labels on the left within each group (optional)
% for g = 1:nGroups
%     for k = 1:nComp
%         text(-0.02*xMax, yPos(g,k)+barH/2, components{k}, ...
%             'HorizontalAlignment','right', ...
%             'VerticalAlignment','middle', ...
%             'FontSize',10);
%     end
% end
% Create dummy handles for legend (so legend shows 3 components)
hLeg = gobjects(nComp,1);
for k = 1:nComp
    hLeg(k) = patch(nan, nan, (1-color_map(k,:))/4.5+color_map(k,:), 'EdgeColor','none');
end
legend(hLeg, components, 'Location','best','Interpreter','tex');

% Make room for the left component labels
% set(gca,'Position',[0.22 0.10 0.73 0.85]); % adjust as needed

%% ====== Draw bracket covering 3 bars + mean value text on the right ======
% Bracket style: a vertical line with short "caps" at top/bottom,
% placed to the right of the longest bar in that group.
bracketPadX = 0.1*xMax; % horizontal offset after max bar in group
capLen      = 0.015*xMax; % cap length
meanPadX    = 0.015*xMax;

for g = 1:nGroups
    vals = dataCell{g};
    xB = max(vals) + bracketPadX;

    yTop = yPos(g,1);                 % first bar y
    yBot = yPos(g,end) + barH;        % last bar bottom edge

    % vertical bracket line
    plot([xB xB],[yTop yBot],'k-','LineWidth',1.2,'HandleVisibility','off');
    % caps
    plot([xB-capLen xB],[yTop yTop],'k-','LineWidth',1.2,'HandleVisibility','off');
    plot([xB-capLen xB],[yBot yBot],'k-','LineWidth',1.2,'HandleVisibility','off');

    % mean text
    meanTxt = sprintf('avg=%.4f', avgVals(g));
    text(xB + meanPadX, (yTop+yBot)/2, meanTxt, ...
         'VerticalAlignment','middle', ...
         'HorizontalAlignment','left', ...
         'FontSize',10, 'Color','k');
end

% ===== Add some top padding in Y =====
yMin = min(yPos(:));
yMax = max(yPos(:)) + barH;

topPad = 0.5;   % 上方留白（单位与 y 相同，可调大/调小）
botPad = 0.5;   % 下方留白（可选）

ylim([yMin - topPad, yMax + botPad]);

%% ====== Cosmetics ======
grid off;
set(gca,'YDir','reverse','Box','on'); % make first group appear at top (like tables)
set(gca,'fontsize',16,'LineWidth',1.5,'TickLength',[0.01, 0.02])

set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');

% title('Relative reconstruction error (component-level)');

% Optional: different colors for x1/x2/x3 within a group
% If you want that, tell me and I'll modify to color-code components.

print(sprintf('./pic/%s.png',filename), '-dpng', '-r400'); % 导出为 PNG 格式
