function createSGMDGui(sig_noisy, fs)
% createSGMDGui 交互式实时更新 SGMD 分解与绘图（带防抖）
% 需要：变量 sig_noisy, fs；路径中含 SGMD 函数
%
% 调用：createSGMDGui(sig_noisy, fs)

    arguments
        sig_noisy {mustBeVector, mustBeNumeric}
        fs (1,1) {mustBePositive}
    end

    % 统一成列向量处理
    x = sig_noisy(:);
    N = numel(x);
    t = (0:N-1).' / fs;

    % 默认参数
    def.threshold_corr = 0.85;
    def.threshold_nmse = 0.01;
    def.tau = 1;   % 固定为 1
    def.debounce = 0.50; % 防抖秒数

    % ---------- UI ----------
    f = uifigure('Name','SGMD 参数探索器（实时）','Position',[100 100 1100 700]);

    gl = uigridlayout(f,[6,6]);
    gl.RowHeight  = {35,35,35,25,'1x','1x'};
    gl.ColumnWidth = {140,120,120,120,120,'1x'};

    % 标题
    lblTitle = uilabel(gl,'Text','SGMD 参数调整与可视化（实时更新）','FontSize',16,'FontWeight','bold');
    lblTitle.Layout.Row = 1;  lblTitle.Layout.Column = [1 6];

    % threshold_corr
    lblCorr = uilabel(gl,'Text','threshold\_corr','HorizontalAlignment','right');
    lblCorr.Layout.Row = 2;  lblCorr.Layout.Column = 1;

    sCorr = uislider(gl,'Limits',[0 1],'Value',def.threshold_corr,'MajorTicks',0:0.1:1);
    sCorr.Layout.Row = 2;  sCorr.Layout.Column = [2 5];

    eCorr = uieditfield(gl,'numeric','Limits',[0 1],'Value',def.threshold_corr,'RoundFractionalValues','off');
    eCorr.Layout.Row = 2;  eCorr.Layout.Column = 6;

    % threshold_nmse
    lblNmse = uilabel(gl,'Text','threshold\_nmse','HorizontalAlignment','right');
    lblNmse.Layout.Row = 3;  lblNmse.Layout.Column = 1;

    sNmse = uislider(gl,'Limits',[0 1],'Value',def.threshold_nmse,'MajorTicks',0:0.1:1);
    sNmse.Layout.Row = 3;  sNmse.Layout.Column = [2 5];

    eNmse = uieditfield(gl,'numeric','Limits',[0 1],'Value',def.threshold_nmse,'RoundFractionalValues','off');
    eNmse.Layout.Row = 3;  eNmse.Layout.Column = 6;

    % 运行按钮 + 信息
    btn = uibutton(gl,'Text','运行 SGMD（可选）','ButtonPushedFcn',@(~,~)runSGMDNow());
    btn.Layout.Row = 4;  btn.Layout.Column = [2 3];

    infoLbl = uilabel(gl,'Text','就绪','FontColor',[0.2 0.2 0.2]);
    infoLbl.Layout.Row = 4;  infoLbl.Layout.Column = [4 6];

    % 上方：原始/重构
    axTop = uiaxes(gl);
    axTop.Layout.Row = 5;  axTop.Layout.Column = [1 6];
    title(axTop,'原始信号与重构信号'); xlabel(axTop,'时间 / s'); ylabel(axTop,'幅值'); grid(axTop,'on');

    % 下方：分量
    axBot = uiaxes(gl);
    axBot.Layout.Row = 6;  axBot.Layout.Column = [1 6];
    title(axBot,'分解分量（最多显示前 8 个）'); xlabel(axBot,'时间 / s'); ylabel(axBot,'幅值'); grid(axBot,'on');

    % 初次绘制原始信号
    cla(axTop); plot(axTop,t,x,'DisplayName','原始'); legend(axTop,'Location','best');

    % ---------- 实时更新的防抖定时器 ----------
    tmr = timer('ExecutionMode','singleShot', ...
                'StartDelay', def.debounce, ...
                'TimerFcn', @(~,~) runSGMDNow());
    % 窗口关闭时清理
    f.CloseRequestFcn = @(src,~) onClose(src, tmr);

    % ---------- 绑定：slider <-> edit 的双向联动 + 实时调度 ----------
    % MATLAB R2019b+ 有 ValueChangingFcn；老版本没有的话删掉两行也行
    if isprop(sCorr,'ValueChangingFcn')
        sCorr.ValueChangingFcn = @(~,e) onSliderChanging(sCorr,e,eCorr);
        sNmse.ValueChangingFcn = @(~,e) onSliderChanging(sNmse,e,eNmse);
    end
    eCorr.ValueChangedFcn  = @(~,~) onEditChanged(sCorr,eCorr);
    eNmse.ValueChangedFcn  = @(~,~) onEditChanged(sNmse,eNmse);

    % 首次自动运行一次
    scheduleRun();

    % ================== 本地函数 ==================
    function onSliderChanging(slider, e, editfld)
        slider.Value = e.Value;
        editfld.Value = e.Value;
        scheduleRun();
    end

    function onEditChanged(slider, editfld)
        slider.Value = editfld.Value;
        scheduleRun();
    end

    function scheduleRun()
        try
            if strcmp(tmr.Running,'on'); stop(tmr); end
            start(tmr);
        catch
            if isempty(tmr) || ~isvalid(tmr)
                tmr = timer('ExecutionMode','singleShot', ...
                            'StartDelay', def.debounce, ...
                            'TimerFcn', @(~,~) runSGMDNow());
            end
            start(tmr);
        end
    end

    function runSGMDNow()
        threshold_corr = sCorr.Value;
        threshold_nmse = sNmse.Value;

        infoLbl.Text = sprintf('运行中：corr=%.3f, nmse=%.3f ...',threshold_corr,threshold_nmse);
        drawnow limitrate;

        try
            % SGCs = SGMD(sig_noisy, fs, def.tau, threshold_corr, threshold_nmse);
            SGCs = SpSSA(sig_noisy, fs, def.tau, 2, ...
                threshold_corr, threshold_nmse, threshold_corr, ...
                'kmeans', true);
            C = normalizeComponents(SGCs, N);
            xhat = sum(C,1).';

            cla(axTop);
            plot(axTop,t,x,'DisplayName','原始'); hold(axTop,'on');
            plot(axTop,t,xhat,'DisplayName','重构');
            hold(axTop,'off');
            xlabel(axTop,'时间 / s'); ylabel(axTop,'幅值'); grid(axTop,'on');
            title(axTop, sprintf('原始 vs 重构（分量数 = %d）', size(C,1)));
            legend(axTop,'Location','best');

            cla(axBot);
            kShow = min(8, size(C,1)-1);
            offs = linspace(0, 0.8, kShow);
            hold(axBot,'on');
            for k = 1:kShow
                plot(axBot, t, C(k,:)' + offs(k), 'DisplayName', sprintf('C%d',k));
            end
            hold(axBot,'off');
            xlabel(axBot,'时间 / s'); ylabel(axBot,'幅值(+offset)'); grid(axBot,'on');
            title(axBot, sprintf('分解分量（显示前 %d 个）', kShow));
            legend(axBot,'Location','eastoutside');

            nmse_val = calcNMSE(x, xhat);
            infoLbl.Text = sprintf('完成：分量=%d，重构NMSE=%.4f（corr=%.3f, nmse=%.3f）',...
                size(C,1), nmse_val, threshold_corr, threshold_nmse);

        catch ME
            infoLbl.Text = ['错误：', ME.message];
            warning(ME.getReport);
        end
    end

end

% ---------- 工具函数 ----------
function C = normalizeComponents(SGCs, N)
    if isnumeric(SGCs)
        [r,c] = size(SGCs);
        if c == N
            C = SGCs;
        elseif r == N
            C = SGCs.';
        else
            error('SGMD 返回的尺寸与信号长度不匹配：[%d x %d] vs N=%d', r, c, N);
        end
    elseif iscell(SGCs)
        K = numel(SGCs);
        C = zeros(K, N);
        for k = 1:K
            ck = SGCs{k}(:);
            if numel(ck) ~= N, error('第 %d 个分量长度不匹配。',k); end
            C(k,:) = ck.';
        end
    elseif isstruct(SGCs)
        fns = fieldnames(SGCs);
        pick = intersect(lower(fns), lower(["SGCs","comps","components","IMFs","modes"]));
        if ~isempty(pick)
            A = SGCs.(pick{1});
            C = normalizeComponents(A, N);
        else
            error('未知的 SGMD 结构体字段，请手动适配。可包含分量矩阵的字段名：SGCs/comps/components/IMFs/modes');
        end
    else
        error('不支持的 SGMD 返回类型：%s', class(SGCs));
    end
end

function v = calcNMSE(x, xhat)
    num = sum((x - xhat).^2);
    den = max(sum(x.^2), eps);
    v = num / den;
end

function onClose(figHandle, tmr)
    try
        if ~isempty(tmr) && isvalid(tmr)
            if strcmp(tmr.Running,'on'); stop(tmr); end
            delete(tmr);
        end
    catch
    end
    delete(figHandle);
end
