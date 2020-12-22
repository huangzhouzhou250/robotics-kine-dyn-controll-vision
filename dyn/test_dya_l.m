% 拉格朗日动力学方程求解测试
%只有 LagrangianDynamics正确
%网址：https://www.jianshu.com/p/6d04539f1cfe
% clc;
% clear;
% syms q1 q2 q3 m1 m2 m3 d2 real
% syms Ix1 Iy1 Iz1 Ixy1 Iyz1 Ixz1 real
% syms Ix2 Iy2 Iz2 Ixy2 Iyz2 Ixz2 real
% syms Ix3 Iy3 Iz3 Ixy3 Iyz3 Ixz3 real
% syms xc1 yc1 zc1 xc2 yc2 zc2 xc3 yc3 zc3 real
% % 
% 
% for i=1:3
%   eval(['Ix',num2str(i),'=','Ix',num2str(i),'+m',num2str(i),'*(yc',num2str(i),'^2+zc',num2str(i),'^2);']);
%   eval(['Iy',num2str(i),'=','Iy',num2str(i),'+m',num2str(i),'*(xc',num2str(i),'^2+zc',num2str(i),'^2);']);
%   eval(['Iz',num2str(i),'=','Iz',num2str(i),'+m',num2str(i),'*(yc',num2str(i),'^2+xc',num2str(i),'^2);']);
%   eval(['Ixy',num2str(i),'=','Ixy',num2str(i),'-m',num2str(i),'*(xc',num2str(i),'*yc',num2str(i),');']);
%   eval(['Ixz',num2str(i),'=','Ixz',num2str(i),'-m',num2str(i),'*(xc',num2str(i),'*zc',num2str(i),');']);
%   eval(['Iyz',num2str(i),'=','Iyz',num2str(i),'-m',num2str(i),'*(yc',num2str(i),'*zc',num2str(i),');']);
% end
% dh_params = [-pi/2, 0,  0, q1; 
%              pi/2, 0, d2, q2;
%              0, 0, 0, q3];
% mass_center = [xc1, yc1, zc1; 
%                xc2, yc2, zc2;
%                xc3, yc3, zc3];
% mass = [m1,m2,m3];
% inertia_1 = [Ix1,  Ixy1, Ixz1;
%              Ixy1, Iy1,  Iyz1;
%              Ixz1, Iyz1, Iz1];
% inertia_2 = [Ix2,  Ixy2, Ixz2;
%              Ixy2, Iy2,  Iyz2;
%              Ixz2, Iyz2, Iz2];
% inertia_3 = [Ix3,  Ixy3, Ixz3;
%              Ixy3, Iy3,  Iyz3;
%              Ixz3, Iyz3, Iz3];
% inertia_tensor(:,:,1) = inertia_1;
% inertia_tensor(:,:,2) = inertia_2;
% inertia_tensor(:,:,3) = inertia_3;
% % 
% [H1,C1,G1] = LagrangianDynamics(dh_params, mass, mass_center, inertia_tensor);



%% puma560测试 
% clc;
% clear;
% syms q1 q2 q3 q4 q5 q6 real
% syms a2 a3 d3 d4 real
% syms m1 m2 m3 m4 m5 m6 real;
% syms Ix1 Iy1 Iz1 Ixy1 Iyz1 Ixz1 real
% syms Ix2 Iy2 Iz2 Ixy2 Iyz2 Ixz2 real
% syms Ix3 Iy3 Iz3 Ixy3 Iyz3 Ixz3 real
% syms Ix4 Iy4 Iz4 Ixy4 Iyz4 Ixz4 real
% syms Ix5 Iy5 Iz5 Ixy5 Iyz5 Ixz5 real
% syms Ix6 Iy6 Iz6 Ixy6 Iyz6 Ixz6 real
% syms xc1 yc1 zc1 xc2 yc2 zc2 xc3 yc3 zc3 real
% syms xc4 yc4 zc4 xc5 yc5 zc5 xc6 yc6 zc6 real
% 
% %alpha a d theta 
% dh_params = [pi/2,  0,  0,  q1;
%                0,  a2,  0,  q2;
%             -pi/2, a3, d3,  q3;
%              pi/2,  0, d4,  q4;
%             -pi/2,  0,  0,  q5;
%                 0,  0,  0,  q6];
% mass_center = [xc1, yc1, zc1;
%                xc2, yc2, zc2;
%                xc3, yc3, zc3;
%                xc4, yc4, zc4;
%                xc5, yc5, zc5;
%                xc6, yc6, zc6;];
% mass = [m1,m2,m3,m4,m5,m6];
% inertia_1 = [Ix1,  Ixy1, Ixz1;
%             Ixy1, Iy1,  Iyz1;
%             Ixz1, Iyz1, Iz1];
% inertia_2 = [Ix2,  Ixy2, Ixz2;
%              Ixy2, Iy2,  Iyz2;
%              Ixz2, Iyz2, Iz2];
% inertia_3 = [Ix3,  Ixy3, Ixz3;
%             Ixy3, Iy3,  Iyz3;
%             Ixz3, Iyz3, Iz3];
% inertia_4 = [Ix4,  Ixy4, Ixz4;
%             Ixy4, Iy4,  Iyz4;
%             Ixz4, Iyz4, Iz4];
% inertia_5 = [Ix5,  Ixy5, Ixz5;
%              Ixy5, Iy5,  Iyz5;
%              Ixz5, Iyz5, Iz5];
% inertia_6 = [Ix6,  Ixy6, Ixz6;
%             Ixy6, Iy6,  Iyz6;
%             Ixz6, Iyz6, Iz6] ;      
% inertia_tensor(:,:,1) = inertia_1;
% inertia_tensor(:,:,2) = inertia_2;
% inertia_tensor(:,:,3) = inertia_3;
% inertia_tensor(:,:,4) = inertia_4;
% inertia_tensor(:,:,5) = inertia_5;
% inertia_tensor(:,:,6) = inertia_6;
% [h,c,g] = LagrangianDynamics(dh_params, mass, mass_center, inertia_tensor);
% mdl_puma560;
% %    eval(['q1(i)=','qr(',num2str(i),')',';']);
% %      q2(i)=eval(['qr(',num2str(i),')',';']);
% %      eval(['qk',num2str(i),'=','qr(i)',';']);
% gc=p560.gravity(3);
% % 
% for i=1:length(qr)
%   link=p560.links(i);
%   eval(['a',num2str(i),'=','link.a',';']);
%   eval(['d',num2str(i),'=','link.d',';']);
% %   eval(['q',num2str(i),'=','qr(i)',';']);
%   eval(['m',num2str(i),'=','link.m',';']);
%   eval(['xc',num2str(i),'=','link.r(1)',';']);
%   eval(['yc',num2str(i),'=','link.r(2)',';']);
%   eval(['zc',num2str(i),'=','link.r(3)',';']);
%   eval(['Ix',num2str(i),'=','link.I(1,1)+link.m*(link.r(2)^2+link.r(3)^2)',';']);
%   eval(['Iy',num2str(i),'=','link.I(2,2)+link.m*(link.r(1)^2+link.r(3)^2)',';']);
%   eval(['Iz',num2str(i),'=','link.I(3,3)+link.m*(link.r(2)^2+link.r(1)^2)',';']);
%   eval(['Ixy',num2str(i),'=','link.I(1,2)+link.m*(link.r(2)*link.r(1))',';']);
%   eval(['Ixz',num2str(i),'=','link.I(1,3)+link.m*(link.r(1)*link.r(3))',';']);
%   eval(['Iyz',num2str(i),'=','link.I(2,3)+link.m*(link.r(2)*link.r(3))',';']);
% end
% dh_params = [pi/2,  0,  0,  q1;
%                0,  a2,  0,  q2;
%             -pi/2, a3, d3,  q3;
%              pi/2,  0, d4,  q4;
%             -pi/2,  0,  0,  q5;
%                 0,  0,  0,  q6];
% mass_center = [xc1, yc1, zc1;
%                xc2, yc2, zc2;
%                xc3, yc3, zc3;
%                xc4, yc4, zc4;
%                xc5, yc5, zc5;
%                xc6, yc6, zc6;];
% mass = [m1,m2,m3,m4,m5,m6];
% inertia_1 = [Ix1,  Ixy1, Ixz1;
%             Ixy1, Iy1,  Iyz1;
%             Ixz1, Iyz1, Iz1];
% inertia_2 = [Ix2,  Ixy2, Ixz2;
%              Ixy2, Iy2,  Iyz2;
%              Ixz2, Iyz2, Iz2];
% inertia_3 = [Ix3,  Ixy3, Ixz3;
%             Ixy3, Iy3,  Iyz3;
%             Ixz3, Iyz3, Iz3];
% inertia_4 = [Ix4,  Ixy4, Ixz4;
%             Ixy4, Iy4,  Iyz4;
%             Ixz4, Iyz4, Iz4];
% inertia_5 = [Ix5,  Ixy5, Ixz5;
%              Ixy5, Iy5,  Iyz5;
%              Ixz5, Iyz5, Iz5];
% inertia_6 = [Ix6,  Ixy6, Ixz6;
%             Ixy6, Iy6,  Iyz6;
%             Ixz6, Iyz6, Iz6] ;      
% inertia_tensor(:,:,1) = inertia_1;
% inertia_tensor(:,:,2) = inertia_2;
% inertia_tensor(:,:,3) = inertia_3;
% inertia_tensor(:,:,4) = inertia_4;
% inertia_tensor(:,:,5) = inertia_5;
% inertia_tensor(:,:,6) = inertia_6;
% 
% [h,c,g] = LagrangianDynamics(dh_params, mass, mass_center, inertia_tensor);
% % h=simplify(h);
% % vpa(h,2)
% for i=1:length(qr)
%     eval(['q',num2str(i),'=','qn(i)',';']);
%     eval(['dq',num2str(i),'=','qn(i)',';']);
% end
% h=subs(h)
% c=subs(c)
% g=subs(g)
% tor=h*qn'+c*qn'+g
%% dyn_lagr_sym验证
% syms q1 q2 q3 m1 m2 m3 d2 real;
% dhtable = [-pi/2, 0,  0, 0; 
%              pi/2, 0, d2, 0;
%              0, 0, 0, 0];
% [D,H,G,fv,fc]=dyn_lagr_sym(dhtable);


% end
% syms a2 a3 d3 d4 real;
% syms q1 q2 q3 q4 q5 q6 real;
% dhtable = [pi/2,  0,  0,  0;
%                0,  a2,  0,  0;
%             -pi/2, a3, d3,  0;
%              pi/2,  0, d4,  0;
%             -pi/2,  0,  0,  0;
%                 0,  0,  0,  0;
% %             
% %             
% [D,H,G,fv,fc]=dyn_lagr_sym(dhtable);
% mdl_puma560;
% %    eval(['q1(i)=','qr(',num2str(i),')',';']);
% %      q2(i)=eval(['qr(',num2str(i),')',';']);
% %      eval(['qk',num2str(i),'=','qr(i)',';']);
% gc=p560.gravity(3);
% 
% for i=1:length(qr)
%   link=p560.links(i);
%   eval(['a',num2str(i),'=','link.a',';']);
%   eval(['d',num2str(i),'=','link.d',';']);
% %   eval(['q',num2str(i),'=','qr(i)',';']);
%   eval(['m',num2str(i),'=','link.m',';']);
%   eval(['xc',num2str(i),'=','link.r(1)',';']);
%   eval(['yc',num2str(i),'=','link.r(2)',';']);
%   eval(['zc',num2str(i),'=','link.r(3)',';']);
%   eval(['Ix',num2str(i),'=','link.I(1,1)',';']);
%   eval(['Iy',num2str(i),'=','link.I(2,2)',';']);
%   eval(['Iz',num2str(i),'=','link.I(3,3)',';']);
%   eval(['Ixy',num2str(i),'=','link.I(1,2)',';']);
%   eval(['Ixz',num2str(i),'=','link.I(1,3)',';']);
%   eval(['Iyz',num2str(i),'=','link.I(2,3)',';']);
% end
% for i=1:length(qn)
%     eval(['q',num2str(i),'=','qn(i)',';']);
%     eval(['dq',num2str(i),'=','qn(i)',';']);
% end
% D=subs(D)
% H=subs(H)
% G=subs(G)
% tor=D*qn'+H*qn'+G
% vpa(tor,6)

%% 动力学算法改进验证
% syms a2 a3 d3 d4 real;
% syms q1 q2 q3 q4 q5 q6 real;
% 
% dhtable= [pi/2,  0,  0,  0;
%                0,  a2,  0,  0;
%             -pi/2, a3, d3,  0;
%              pi/2,  0, d4,  0;
%             -pi/2,  0,  0,  0;
%                 0,  0,  0,  0];
%             
%             
% [D,H,G,fv,fc]=dyn_lagr_sym(dhtable);

%% 两连杆动力学算法
clc
clear
syms l1 l2 gc real;
syms q1 q2 dq1 dq2 ddq1 ddq2 real;
g=[gc,0,0,0]';
dhtable=[0 l1 0 0;
    0 l2 0,0;];
[D,H,G,fv,fc]=dyn_lagr_sym(dhtable,g);

Ix1=0;Iy1=0;Iz1=0;Ixx1=0;Iyy1=0;Izz1=0;
Ix2=0;Iy2=0;Iz2=0;Ixx2=0;Iyy2=0;Izz2=0;

% l1=0.25;l2=0.2;
xc1=l1;xc2=l2;
yc1=0;yc2=0;
zc1=0;zc2=0;
% 
% m1=2;m2=1;
% 
% q1=pi/2;q2=pi/2;
% dq1=0.1;dq2=0.2;
% ddq1=1;ddq2=2;
% gc=9.8;
D=subs(D)
H=subs(H)
G=subs(G)

tor=D*[ddq1;ddq2]+H*[dq1;dq2]+G;

% l1=0.25;l2=0.2;
% L1= Link('revolute', 'd', 0,     'a',l1,        'alpha',0,'offset',0,'r',[l1 0 0]','m',2);
% L2= Link('revolute', 'd', 0,         'a',l2,    'alpha', 0,'offset',0,'r',[l2 0 0],'m',1);
% tlink=SerialLink([L1,L2],'gravity',[-gc 0 0]);
% q=[q1,q2];dq=[dq1,dq2];ddq=[ddq1,ddq2];
% tlink.rne(q,dq,ddq)
% vpa(tor,4)



%结果
D = [ 4*l1^2*m1 + l1^2*m2 + 4*l2^2*m2 + 4*l1*l2*m2*cos(q2), 4*m2*l2^2 + 2*l1*m2*cos(q2)*l2;
                       4*m2*l2^2 + 2*l1*m2*cos(q2)*l2,                      4*l2^2*m2];
 H =[ -2*dq2*l1*l2*m2*sin(q2), -2*l1*l2*m2*sin(q2)*(dq1 + dq2);
 2*dq1*l1*l2*m2*sin(q2),                               0]  ;
G=[ gc*m2*(l2*sin(q1 + q2) + l1*sin(q1)) + gc*l2*m2*sin(q1 + q2) + 2*gc*l1*m1*sin(q1);
                                                           2*gc*l2*m2*sin(q1 + q2)];

