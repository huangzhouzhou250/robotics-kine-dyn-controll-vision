%基于指数积公式的机器人运动学仿真平台

%% 绘制基础框架
%绘制框架fig
set(0,'Units','pixels')
dim = get(0,'ScreenSize');
fig_1 = figure('doublebuffer','on','Position',[0,35,dim(3)-150,dim(4)-100],...
    'MenuBar','none','Name','POE Robot Drawing',...
    'NumberTitle','off');
hold on;

light    %灯光                           
daspect([1 1 1])                   
view(135,25)   %视角
xlabel('X'),ylabel('Y'),zlabel('Z');
title('指数积机器人仿真');
axis([-800 800 -800 800 -100 1100]);
%绘制图形边框 
plot3([-800,800],[-800,-800],[-100,-100],'k')
plot3([-800,-800],[-800,800],[-100,-100],'k')
plot3([-800,-800],[-800,-800],[-100,1100],'k')
plot3([-800,800],[-800,-800],[1100,1100],'k')
plot3([-800,-800],[-800,800],[1100,1100],'k')
plot3([-800,-800],[800,800],[-100,1100],'k')

%绘制主要功能区 
kine = uipanel(fig_1,...
    'Position',[0.05 0.15 0.15 0.35],...
    'Title','运动学','FontSize',11);

