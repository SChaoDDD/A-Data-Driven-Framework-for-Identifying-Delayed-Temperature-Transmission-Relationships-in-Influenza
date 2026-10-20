%% SDtemp 温度分析
% 需求：从 2013/9/30 开始往后 2200 天，计算
%   1) 当日温度 T(t)
%   2) 加权温度：当日及前 7 天（共 8 天）加权，权重 w(1)~w(8) 对应 t, t-1, ..., t-7
%   3) 平均温度变化：(T(t) - T(t-30)) / 30
% 输出 CSV：时间, 温度, 加权温度, 平均变化温度

clear; clc; close all;

%% 0. 参数设置
inputFile  = 'C:\Users\XLYC\Desktop\流感代码修改\AAAinfluenzacode\matlab\IA';   % 输入数据路径（按需修改）
outputFile = 'C:\Users\XLYC\Desktop\流感代码修改\AAAinfluenzacode\matlab\IA\IAtemp_result.csv'; % 输出结果路径

startDate  = datetime(2013, 9, 30);      % 起始日
nDays      = 2200;                        % 往后天数
L_mean     = 30;                          % 平均温度变化窗口
L_wgt      = 7;                           % 加权窗口（前 7 天，含当日共 8 天）

% 权重：w(1)=今日, w(2)=前1天, ..., w(8)=前7天
w = [0.5629, 0.1257, 0.1178, 0.0906, 0.0512, 0.0258, 0.0147, 0.0113].';

%% 1. 读取数据
% 文件格式：date,tmean_avg_c （表头行）
T = readtable('IAtemp.csv');
T.Date = datetime(T.date, 'InputFormat', 'yyyy/M/d');  % 解析日期
T.date = [];                                            % 删除原字符列
T = sortrows(T, 'Date');                                % 按日期升序

dates_all = T.Date;
temp_all  = T.MEAN;

%% 2. 构建等间距的日序列并填充缺失
dVec = (dates_all - dates_all(1));
dVec = days(dVec);
if isnat(dates_all(1))
    error('首行日期解析失败，请检查日期格式。');
end
% 以首日为准建立整数日索引
allDates = (dates_all(1):days(1):dates_all(end)).';
nAll     = numel(allDates);
tempMap  = NaN(nAll, 1);
[~, idxInAll] = ismember(dates_all, allDates);
tempMap(idxInAll) = temp_all;

% 缺失日用前后有效值的线性插值补齐（保证 t-30、t-7 均有值）
if any(isnan(tempMap))
    valid = find(~isnan(tempMap));
    tempMap = interp1(valid, tempMap(valid), 1:nAll, 'linear', 'extrap');
end

%% 3. 定位起始索引
[~, startIdx] = ismember(startDate, allDates);
if startIdx == 0
    error('起始日期 %s 不在数据范围内。', datestr(startDate));
end
% 需要 t-30 到 t 共 31 个值，故起始前至少要有 30 天数据
if startIdx <= L_mean
    error('起始日期前不足 %d 天数据，无法计算。', L_mean);
end

endIdx = startIdx + nDays - 1;
if endIdx > nAll
    error('数据只到 %s，不足 %d 天。', datestr(allDates(end)), nDays);
end

%% 4. 主计算循环
selDates = allDates(startIdx:endIdx);
n = numel(selDates);

T_series  = zeros(n, 1);   % 当日温度
T_weight  = zeros(n, 1);   % 加权温度
T_change  = zeros(n, 1);   % 平均温度变化

for k = 1:n
    cur = startIdx + k - 1;       % 当前在 tempMap 中的位置
    T_series(k) = tempMap(cur);

    % --- 加权温度：cur, cur-1, ..., cur-L_wgt (共 L_wgt+1 = 8 个点)
    seg = tempMap(cur - (0:L_wgt)).';   % 长度 8，[今日; 前1天; ...; 前7天]
    T_weight(k) = seg*w;           % 加权求和

    % --- 平均温度变化：(T(t) - T(t-30)) / 30
    T_change(k) = (tempMap(cur) - tempMap(cur - L_mean)) / L_mean;
end

%% 5. 组装并输出 CSV
outTable = table(...
    selDates, ...
    round(T_series, 4), ...
    round(T_weight, 4), ...
    round(T_change, 4), ...
    'VariableNames', {'Time', 'Temperature', 'WeightedTemperature', 'MeanChangeTemperature'});

writetable(outTable,'result.csv');

fprintf('计算完成，共 %d 行，已写入：\n%s\n', n, outputFile);
fprintf('时间范围：%s ~ %s\n', datestr(selDates(1)), datestr(selDates(end)));
