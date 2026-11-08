function idxB = alignSignalsByCorr(A, B)
% alignSignalsByCorr  根据列之间的相关系数，将 B 的列重排以对齐 A
%
%   idxB = alignSignalsByCorr(A, B)
%
%   输入：
%       A, B : 大小为 N x M 的矩阵，每列为一个信号。
%               要求 A 和 B 行数、列数都相同。
%
%   输出：
%       idxB : 1 x M 的索引向量，使得
%               B_aligned = B(:, idxB);
%              之后 A(:, i) 与 B_aligned(:, i) 是匹配信号。
%
%   说明：
%       - 使用皮尔逊相关系数衡量信号相似度；
%       - 使用简单贪心算法：每次选取当前相关系数最大的未匹配列对；
%       - 在大多数噪声不太大的场景下效果够用，如果需要全局最优可换成匈牙利算法。

    % --- 基本检查 ---
    [nA, mA] = size(A);
    [nB, mB] = size(B);

    if nA ~= nB
        error('A 和 B 的行数（信号长度）必须相同。');
    end
    if mA ~= mB
        error('当前实现要求 A 和 B 的列数（信号个数）相同。');
    end

    M = mA;

    % --- 计算相关系数矩阵 S（M x M）---
    S = zeros(M, M);
    for i = 1:M
        ai = A(:, i);
        for j = 1:M
            bj = B(:, j);
            c = corr(ai, bj);   % 皮尔逊相关系数
            if isnan(c)
                c = 0;          % 防止 NaN 破坏后续 max
            end
            S(i, j) = c;
        end
    end

    % --- 贪心匹配：每次选择当前 S 中最大值 ---
    S_work = S;                 % 工作副本
    idxB   = zeros(1, M);       % 最终索引：A 的第 i 列匹配 B 的第 idxB(i) 列

    for k = 1:M
        % 每行最大值及其列索引
        [maxRowVals, bestColIdx] = max(S_work, [], 2);
        % 在所有行中再取最大的那个
        [maxVal, iA] = max(maxRowVals);

        if ~isfinite(maxVal)
            % 没有可用的匹配了（理论上不会出现，除非 M=0）
            break;
        end

        jB = bestColIdx(iA);  % A 的第 iA 列匹配 B 的第 jB 列

        idxB(iA) = jB;

        % 屏蔽掉已经匹配过的行、列
        S_work(iA, :) = -Inf;
        S_work(:, jB) = -Inf;
    end
end
