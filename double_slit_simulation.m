% 双缝夫琅禾费衍射/干涉光学仿真
% 运行本脚本即可生成二维光强图和中心截面曲线
clear; close all; clc;

% ---- 物理参数（SI单位） ----
lambda = 633e-9;      % 波长：633 nm（红光）
a = 0.08e-3;          % 单缝宽度：0.08 mm
d = 0.35e-3;          % 双缝中心距：0.35 mm
z = 1.2;              % 传播距离：1.2 m
I0 = 1;               % 入射归一化强度

% ---- 观察屏网格 ----
N = 900;
x = linspace(-0.012, 0.012, N);   % 水平范围 ±12 mm
y = linspace(-0.008, 0.008, round(N*2/3));
[X,Y] = meshgrid(x,y);

% 夫琅禾费双缝强度：单缝包络 × 双缝干涉项
thetaX = atan2(X,z);
thetaY = atan2(Y,z);
beta = pi*a*sin(thetaX)/lambda;
alpha = pi*d*sin(thetaX)/lambda;
sincTerm = ones(size(beta));
idx = abs(beta) > 1e-12;
sincTerm(idx) = sin(beta(idx))./beta(idx);
I = I0 .* sincTerm.^2 .* cos(alpha).^2;

% ---- 绘图与保存 ----
fig = figure('Color','w','Position',[100 100 1100 700]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
nexttile;
imagesc(x*1e3,y*1e3,I); axis image; set(gca,'YDir','normal');
xlabel('屏幕 x (mm)'); ylabel('屏幕 y (mm)');
title('双缝夫琅禾费衍射/干涉光强'); colormap hot; colorbar;
nexttile;
plot(x*1e3,I(round(end/2),:),'r','LineWidth',1.4); grid on;
xlabel('屏幕 x (mm)'); ylabel('归一化光强');
title('中心截面光强曲线'); xlim([min(x) max(x)]*1e3);

sgtitle(sprintf('\\lambda = %.0f nm, a = %.2f mm, d = %.2f mm, z = %.2f m', ...
    lambda*1e9,a*1e3,d*1e3,z));
outFile = fullfile(pwd,'double_slit_result.png');
exportgraphics(fig,outFile,'Resolution',150);
save(fullfile(pwd,'double_slit_result.mat'),'lambda','a','d','z','x','y','I');
fprintf('仿真完成！\n最大归一化光强: %.3f\n结果图: %s\n',max(I,[],'all'),outFile);
