function [fig,ax,dcm_obj] = fig_open(width, height)
if nargin < 2
    width = 600;
    height = 500;
end
set(0,'defaultfigurecolor',[1 1 1]) % white background
set(0,'defaultaxesfontname','Arial') % beautify the axes a bit
scrn_size = get(0,'ScreenSize'); % get size of screen
shrink_pct = 0.1; % shrink the figure by 10%
%
fig = figure('Visible','on','DefaultAxesFontSize',6,'Position',...
    [scrn_size(1)+(scrn_size(3)*shrink_pct) scrn_size(2)+(scrn_size(4)*shrink_pct)...
    scrn_size(3)-(scrn_size(3)*2*shrink_pct) scrn_size(4)-(scrn_size(4)*2*shrink_pct)]); % shrinking the figure
% fig = figure('Visible','on','DefaultAxesFontSize',20);
%
ax = gca;
dcm_obj = datacursormode(fig); % enable control of datatips
% data tips
set(ax,'fontsize',3,'LineWidth',1.5,'TickLength',[0.025, 0.02])
% set(ax,'Color',[0.8 0.8 0.8],'gridcolor',[1 1 1],'gridalpha',0.9) % set the axis color
% to compliment the white background
%
hold all
% grid on % grid for more contrast in the axes
% set(gca, 'XGrid', 'off');   % 关闭竖向网格线
% set(gcf, 'Position', [0, 0, width, height]);
end