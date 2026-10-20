clc;
clear;
close all;

% T 的取值范围
T = linspace(-27, 31, 100);

% 固定的不同 dT 值
dT_values = [-10/30, -5/30, 0/30, 5/30, 10/30];

% 创建图像
figure;
hold on;
grid on;
box on;

% 对每个固定 dT 计算函数值
for i = 1:length(dT_values)
    
    dT = dT_values(i);
    
    y = exp(1.0408957382.*(0.0000057972.*T.^3-0.0011048007.*T.^2+0.0050645368.*T-0.4386011040.*dT+0.1641361835))./7;
    
    plot(T, y, 'LineWidth', 2, ...
        'DisplayName', sprintf('dT = %.1f', dT));
end

% 坐标轴设置
xlabel('T');
ylabel('f(T, dT)');


% 显示图例
legend('Location', 'best');

% 字体设置
set(gca, 'FontSize', 12);