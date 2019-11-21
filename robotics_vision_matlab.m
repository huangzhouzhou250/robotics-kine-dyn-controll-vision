%% 机器人视觉控制仿真学习的代码内容
%% 重要代码

%  rtbdemo 
% 运行 startup_rvc 
% SerialLink.plot3d 机器人模型的三维演示
%安装视觉工具包后，运行动力学仿真程序会出问题
%%
% TB = SE3(1, 2, 0) * SE3.Rz(pi/2);
% vb = [0.2 0.3 0 0 0 0.5]';
% va = TB.velxform * vb;

%3.3 Creating Time-Varying Pose

% tpoly(0, 1, 50);
% [s,sd,sdd] = tpoly(0, 1, 50);
% tpoly(0, 1, 50, 0.5, 0);
% mean(sd) / max(sd)  %角速度利用率不够

% lspb(0, 1, 50);
% [s,sd,sdd] = lspb(0, 1, 50);
% lspb(0, 1, 50, 0.025);
% lspb(0, 1, 50, 0.035);


% q = mtraj(@lspb, [0 2], [1 -1], 50);
% plot(q)

% via = SO2(30, 'deg') * [-1 1; 1 1; 1 -1; -1 -1]';
% q0 = mstraj(via(:,[2 3 4 1])', [2,1], [], via(:,1)', 0.2, 0);
% plot(q0(:,1), q0(:,2))

% q2 = mstraj(via(:,[2 3 4 1])', [2,1], [], via(:,1)', 0.2, 2);

% R0 = SO3.Rz(-1) * SO3.Ry(-1);
% R1 = SO3.Rz(1) * SO3.Ry(1);
% rpy0 = R0.torpy(); rpy1 = R1.torpy();
% rpy = mtraj(@tpoly, rpy0, rpy1, 50);
% SO3.rpy( rpy ). animate;
% 
% q0 = R0. UnitQuaternion; q1 = R1.UnitQuaternion;
% q = interp(q0, q1, 50);
% q.animate

% T0 = SE3([0.4, 0.2, 0]) * SE3.rpy(0, 0, 3);
% T1 = SE3([-0.4, -0.2, 0.3]) * SE3.rpy(-pi/4, pi/4, -pi/2);
% interp(T0, T1, 0.5)
% Ts = ctraj(T0, T1, 50);

%3.4 · Application: Inertial Navigation
% ex_tumble;
% attitude(1) = UnitQuaternion();
% for k=1:numcols(w)-1
% attitude(k+1) = attitude(k) .* UnitQuaternion.omega( w(:,k)*dt );
% end
% attitude.animate('time', t);
% mplot(t, attitude.torpy() )


%% III Arm-Type Robots
%%chapter 7 Robot Arm Kinematics
% import ETS2.*
% a1 = 1;
% import ETS3.*
% L1 = 0; L2 = -0.2337; L3 = 0.4318; L4 = 0.0203; L5 = 0.0837; L6 = 0.4318;


% Denavit-Hartenberg Parameters
%  L = Revolute('a', 1);
%  L.A(0.5);
%  L.type;
%  L.A(0);
%  robot = SerialLink( [ Revolute('a', 1) Revolute('a', 1) ],'name', 'my robot');
%  robot.fkine([30 40], 'deg');
% robot.edit

%  models;%b包含的机器人模型
% %  ABB, IRB140, 6DOF, modified_DH (mdl_irb140_mdh)
% % ABB, IRB140, 6DOF, standard_DH (mdl_irb140)
% % ABB, S4_2.8, S4 2.8m reach version, 6DOF, standard_DH (mdl_S4ABB2p8)
% % Adept, Cobra600, 4DOF, standard_DH (mdl_cobra600)
% % Aldebaran, NAO, humanoid, 4DOF, standard_DH (mdl_nao)
% % Baxter, Rethink Robotics, 7DOF, standard_DH (mdl_baxter)
% % Fanuc, AM120iB/10L, 6DOF, standard_DH (mdl_fanuc10L)
% % Fanuc, M16, 6DOF, standard_DH (mdl_M16)
% % Kinova, Jaco, 6DOF, standard_DH (mdl_jaco)
% % Kinova, Mico, 6DOF, standard_DH (mdl_mico)
% % Kuka, KR5, 6DOF, standard_DH (mdl_KR5)
% % Kuka, LWR, 7DOF, standard_DH (mdl_LWR)
% % Motoman, HP6, 6DOF, standard_DH (mdl_motomanHP6)
% % Stanford, Stanford Arm, prismatic, 6DOF, standard_DH (mdl_stanford)
% % Stanford, Stanford arm, prismatic, 6DOF, modified_DH (mdl_stanford_mdh)
% % Trossen Robotics, PhantomX Pincher, 4DOF, standard_DH (mdl_phantomx)
% % Unimation, Puma560, dynamics, 6DOF, modified_DH (mdl_puma560akb)
% % Unimation, Puma560, dynamics, 6DOF, standard_DH (mdl_puma560)
% % Unimation, Puma560, on XY base, redundant, 8DOF, standard_DH (mdl_p8)
% % Universal Robotics, UR10, 6DOF, standard_DH (mdl_ur10)
% % Universal Robotics, UR3, 6DOF, standard_DH (mdl_ur3)
% % Universal Robotics, UR5, 6DOF, standard_DH (mdl_ur5)
% % generic, 6DOF, standard_DH (mdl_offset6)
% % generic, 6DOF, standard_DH (mdl_simple6)
% % generic, ball shape, hyper redundant, 50DOF, standard_DH (mdl_ball)
% % generic, coil, hyper redundant, 50DOF, standard_DH (mdl_coil)
% % generic, planar, 1DOF, standard_DH (mdl_onelink)
% % generic, planar, 1DOF, standard_DH (mdl_planar1)
% % generic, planar, 2DOF, modified_DH (mdl_twolink_mdh)
% % generic, planar, 2DOF, standard_DH (mdl_planar2)
% % generic, planar, 2DOF, symbolic, standard_DH (mdl_planar2_sym)
% % generic, planar, 3DOF, standard_DH (mdl_planar3)
% % generic, planar, dynamics, 2DOF, standard_DH (mdl_twolink)
% % generic, planar, dynamics, 2DOF, symbolic, standard_DH (mdl_twolink_sym)
% % hyper redundant, 10DOF, standard_DH (mdl_hyper3d)
% % planar, hyper redundant, 10DOF, standard_DH (mdl_hyper2d)

% mdl_irb140

% a1 = 1; a2 = 1;
% TE0 = SE2(a1+a2, 0, 0);
% S1 = Twist( 'R', [0 0] );
% S2 = Twist( 'R', [a1 0] );
% TE = S1.T(30*pi/180) * S2.T(40*pi/180) * TE0.T

% 7.1.2.3 l 6-Axis Industrial Robot
% mdl_puma560;
% p560;
% p560.plot(qz);
% TE = p560.fkine(qz);
% p560.tool = SE3(0, 0, 0.2);
% p560.fkine(qz)
% p560.base = SE3(0, 0, 30*0.0254);
% q=jtraj(qz,qr,10);
% T = p560.fkine(q);

% 7.2lInverse Kinematics
% mdl_puma560;
% qn;
% T = p560.fkine(qn);
% qi = p560.ikine6s(T);
% p560.ikine6s( SE3(3, 0, 0) );%封闭解
% q = [0 pi/4 pi 0.1 0 0.2];
% p560.ikine6s(p560.fkine(q), 'ru');
% q(4)+q(6);%如何在这设置的角度范围
% T = p560.fkine(qn);
% qi = p560.ikine(T);
% p560.fkine(qi);
% qi = p560.ikine(T, 'q0', [0 0 3 0 0 0]);%

%Under-Actuated Manipulator
% mdl_cobra600;
% c600;
% T = SE3(0.4, -0.3, 0.2) * SE3.rpy(30, 40, 160, 'deg');
%  q = c600.ikine(T, 'mask', [1 1 1 0 0 1]) 
% Ta = c600.fkine(q);‘

%7.2.2.4l Redundant Manipulator
%  mdl_baxter;
%  left;
%  q = left.ikine(TE);
% TE = SE3(0.8, 0.2, -0.2) * SE3.Ry(pi);
% left.fkine(q).print('xyz')
% left.plot(q)

%7.3 Trajectories
% mdl_puma560
% T1 = SE3(0.4, 0.2, 0) * SE3.Rx(pi);
% T2 = SE3(0.4, -0.2, 0) * SE3.Rx(pi/2);
% q1 = p560.ikine6s(T1);
% q2 = p560.ikine6s(T2);
% t = [0:0.05:2]';
% q = jtraj(q1, q2, t);
% [q,qd,qdd] = jtraj(q1, q2, t);
% q = p560.jtraj(T1, T2, t);%多种求解关节插值的方法
% p560.plot(q);
% qplot(t, q);
% T = p560.fkine(q);
% p = T.transl;
% 
% % 7.3.2Cartesian Motion p214
% Ts = ctraj(T1, T2, length(t));
% plot(t, Ts.transl);
% qc = p560.ikine6s(Ts);

% sl_jspace
% T1 = SE3(0.5, 0.3, 0.44) * SE3.Ry(pi/2);
% T2 = SE3(0.5, -0.3, 0.44) * SE3.Ry(pi/2);
%  t = [0:0.05:2]';
% Ts = ctraj(T1, T2, length(t));
% qc = p560.ikine6s(Ts);
% 
% %7.3.5l Configuration Change
% T = SE3(0.4, 0.2, 0) * SE3.Rx(pi);
% qr = p560.ikine6s(T, 'ru');
% ql = p560.ikine6s(T, 'lu');
% q = jtraj(qr, ql, t);
% p560.plot(q);

%%7.4lAdvanced Topics

%%7.5l Applications
% load hershey
% B = hershey{'B'}
% B.stroke  %获取字母的轨迹
% path = [ 0.25*B.stroke; zeros(1,numcols(B.stroke))]; %调整合适大小
% k = find(isnan(path(1,:))); % 寻找其中的过渡点，在原轨迹中，过渡点为NAN
% path(:,k) = path(:,k-1); path(3,k) = 0.2;%使过渡点为悬空暂留
% traj = mstraj(path(:,2:end)', [0.5 0.5 0.5], [], path(:,1)',0.02, 0.2);%线性与多次线混合
% % traj = mstraj(path(:,2:end)'{途径点}, [0.5 0.5 0.5]{速度限制}, []{每段时间}, path(:,1)'{起始位置},0.02{时间步长}, 0.2{加速度转换});
% numrows(traj) * 0.02;
% plot3(traj(:,1), traj(:,2), traj(:,3));
% Tp = SE3(0.6, 0, 0) * SE3(traj) * SE3.oa( [0 1 0], [0 0 -1]);
% q = p560.ikine6s(Tp);
% p560.plot(q)

% A Simple Walking Robot
% s = 'Rz(q1) Rx(q2) Ty(L1) Rx(q3) Tz(L2)';
% dh = DHFactor(s);%运行 startup_rvc 才能使用
% dh.tool;
% dh.command('leg');
% L1 = 0.1; L2 = 0.1;
% leg = eval( dh.command('leg') );
% transl( leg.fkine([0,0,0]) );
% leg.plot([0,0,0], 'nobase', 'noshadow', 'notiles');
% set(gca, 'Zdir', 'reverse'); view(137,48);
% transl( leg.fkine([0.2,0,0]) );
% 
% xf = 50; xb = -xf; y = 50; zu = 20; zd = 50;
% path = [xf y zd; xb y zd; xb y zu; xf y zu; xf y zd] * 1e-3;
% qcycle = leg.ikine( SE3(p), 'mask', [1 1 1 0 0 0] );
% leg.plot(qcycle, 'loop');%loop用于循环画图
% 
% W = 0.1; L = 0.2;
% legs(1) = SerialLink(leg, 'name', 'leg1');
% legs(2) = SerialLink(leg, 'name', 'leg2', 'base', SE3(-L, 0, 0));
% legs(3) = SerialLink(leg, 'name', 'leg3', 'base', SE3(-L, -W, 0)*SE3.Rz(pi));
% legs(4) = SerialLink(leg, 'name', 'leg4', 'base', SE3(0, -W, 0)*
% SE3.Rz(pi));%确定四个腿的位置
% clf; k = 1;
% while 1
% legs(1).plot( gait(qcycle, k, 0, false) );
% if k == 1, hold on; end
% legs(2).plot( gait(qcycle, k, 100, false) );
% legs(3).plot( gait(qcycle, k, 200, true) );
% legs(4).plot( gait(qcycle, k, 300, true) );
% drawnow
% k = k+1;
% end

%% 8 Manipulator Velocity

% mdl_planar2_sym
% p2
% syms q1 q2 real
% TE = p2.fkine( [q1 q2] );
% p = TE.t;p = p(1:2)
% J = jacobian(p, [q1 q2]);
% 
% 
% p560.jacobe(qn)
% J = p560.jacob0(qn)
% 
% rpy2jac(0.1, 0.2, 0.3);
% p560.jacob0(qn, 'eul');
% 
% % 8.2 l Jacobian Condition and Manipulability
% J = p560.jacob0(qr)
% rank(J) %判断雅克比矩阵秩
% jsingu(J) %判断雅克比矩阵降秩的原因
% 
% qns = qr; qns(5) = 5 * pi/180
% J=p560.jacob0(qns);
% qd = inv(J)*[0 0 0.1 0 0 0]' ;
% det(J)
% cond(J) %条件数
% qd = inv(J)*[0 0 0 0 0.2 0]';
% 
% %8.2.2l Manipulability
% mdl_planar2
% p2.vellipse([30 40], 'deg')
% p2.teach('callback', @(r,q) r.vellipse(q), 'view', 'top')
% 
% J = p560.jacob0(qns);
% J = J(1:3, :);
% plot_ellipse(J*J')
% p560.vellipse(qns, 'rot')
% p560.vellipse(qns, 'trans');

% 8.3 Resolved-Rate Motion Control
% sl_rrmc
% r = sim( 'sl_rrmc');
% t = r.find('tout');%获取仿真相关数据
% q = r.find('yout');
% T = p560.fkine(q);
%  xyz = transl(T);
% mplot(t, xyz(:,1:3))
% 
%  mplot(t, q(:,1:3))
%  
%  sl_rrmc2 %仔细研究仿真的内容
%  
% %  8.3.l Jacobian Singularity
% mdl_planar2
% qn = [1 1];
% 
% J = p2.jacob0(qn)
% qd = pinv(J) * [0.1 0 0 0 0 0]'
% xd = J*qd;
% Jxy = J(1:2,:);
% qd = inv(Jxy)* [0.1 0]';
% 
% 
% mdl_baxter
% TE = SE3(0.8, 0.2, -0.2) * SE3.Ry(pi);
% q = left.ikine(TE)
% J= jacob0(left, q);
% xd = [0.2 0.2 0.2 0 0 0]';
% qd = pinv(J) * xd;
% qd'
% N = null(J)
% norm( J * N(:,1))
% qd_null = [0 0 0 0 1 0 0]';

% 8.5 Force Relationships
%  tau = p560.jacob0(qn)' * [20 0 0 0 0 0]'; 
%  p2.fellipse([30 40], 'deg')
%  p2.teach(qn, 'callback', @(r,q) r.fellipse(q), 'view', 'top')   %help teach查看说明
 %可以参考来做界面的绘制
 
%  8.6  Inverse Kinematics: a General Numerical Approach
% ikine %求解方法

% 8.7 Advanced Topics
%  p560.teach(qn,'workspace',[-1 1 -1 1 -1 1]) %设置工作空间

%% 9 Dynamics and Control
%  9.1 Independent Joint  Control 
% mdl_twolink_sym
% syms q1 q2 q1d q2d q1dd q2dd real
% tau = twolink.rne([q1 q2], [q1d q2d], [q1dd q2dd]);
% 
% tf = p560.jointdynamics(qn);
% 
% vloop_test %速度环控制策略 研究一下
% sim('vloop_test')
% 
%  ploop_test %速度环控制，缺少G（电机减速比）%需要研究一下
 
%  9.2 Rigid-Body Equations of Motion
% mdl_puma560
% Q = p560.rne(qn, qz, qz);    %动力学求解，牛顿欧拉法，研究一下
% 
% q = jtraj(qz, qr, 10);
% Q = p560.rne(q, 0*q, 0*q);
% 
% p560.rne(qn, [1 0 0 0 0 0], qz, 'gravity', [0 0 0]);
% p560.links(1).dyn %连杆的动力学参数
% 
% gravload = p560.gravload(qn);
% p560.gravity' %可以通过赋值来改变
% p560.gravity = p560.gravity/6;
% p560.base = SE3.Rx(pi);
% p560.gravload(qn)
% 
% Q = p560.gravload(qs) %第二关节的受力最大
%  [Q2,Q3] = meshgrid(-pi:0.1:pi, -pi:0.1:pi);
%  for i=1:numcols(Q2),
%  for j=1:numcols(Q3);
%  g = p560.gravload([0 Q2(i,j) Q3(i,j) 0 0 0]);
%  g2(i,j) = g(2);
%  g3(i,j) = g(3);
%  end
%  end
%  surfl (Q2, Q3, g2); surfl (Q2, Q3, g3);

% M = p560.inertia(qn);
% C = p560.coriolis(qn, qd)
%  [Q2,Q3] = meshgrid(-pi:0.1:pi, -pi:0.1:pi);
%  for i=1:numcols(Q2)
%      for j=1:numcols(Q3)
%          M = p560.inertia([0 Q2(i,j) Q3(i,j) 0 0 0]);
%          M11(i,j) = M(1,1);
%          M22(i,j) = M(2,2);
%      end
% end
% surfl (Q2, Q3, M11); 
% surfl (Q2, Q3, M22);
% max(M11(:)) / min(M11(:))
% 
% qd = [0 0 1 0 0 0];
% 
% C*qd'
% 
% p560.payload(2.5, [0 0 0.1]);
% M_loaded = p560.inertia(qn);
% M_loaded ./ M
% p560.gravload(qn) ./ gravload
% 
% [Q,Wb] = p560.rne(qn, qz, qz);
% Wb'
% sum([p560.links.m])
% 
% %力椭球
% J = p560.jacob0(qn);
% M = p560.inertia(qn);
% Mx = (J * inv(M) * inv(M)' * J');
% Mx = Mx(1:3, 1:3);
% sqrt(eig(Mx))

% 9.3 Forward Dynamics
%qdd = p560.accel(q, qd, Q) %计算加速度
% sl_ztorque   %动力学仿真 无外力时仿真 有意思
% r = sim('sl_ztorque');
% t = r.find('tout');
% q = r.find('yout');
% plot(t, q(:,1:3))

% 9.4 Rigid-Body Dynamics Compensation  （研究一下两个力控仿真）
% mdl_puma560
%  sl_fforward
% r = sim('sl_fforward'); %前馈力控
% 
% sl_ctorque  %反馈力控
% r = sim( 'sl_ctorque');
% 
% sl_opspace %操作空间力控 

% 9.5 Applications
% sl_sea

%% Part IV Computer Vision
% 10 Light and Color
% 10.1 Spectral Representation of Light
% lambda = [300:10:1000]*1e-9;
% for T=3000:1000:6000
%  plot( lambda, blackbody(lambda, T)); 
%  hold all
% end
%  10.3 Advanced Topics
% lambda = [400:10:700]*1e-9;

%% 11 Image Formation
% 11.1 Perspective Camera
% 
% cam = CentralCamera('focal', 0.015);
% P = [0.3, 0.4, 3.0]';
% cam.project(P)
% cam.project(P, 'pose', SE3(-0.5, 0, 0) )
% 
% cam = CentralCamera('focal', 0.015, 'pixel', 10e-6,...
% 'resolution', [1280 1024], 'centre', [640 512], 'name', 'mycamera')
% cam.project(P) %依据上述参数求解P在像平面中的坐标
% 
% cam.K %相机内部变换矩阵，涉及内部参数
% 
% cam.C %相机的参数矩阵 C
% 
% cam.fov() * 180/pi %镜头观看的角度
% 
% P = mkgrid(3, 0.2, 'pose', SE3(0, 0, 1.0)); %构造栅格数据
% cam.project(P)
% cam.plot(P)
% Tcam = SE3(-1,0,0.5)*SE3.Ry(0.9);
% cam.plot(P, 'pose', Tcam) %以特定的视角观察这些点
% 
% cam.project([1 0 0 0]', 'pose', Tcam)
% p = cam.plot(P, 'pose', Tcam)
% 
% cube = mkcube(0.2, 'pose', SE3(0, 0, 1) ); %绘制正方体
% cam.plot(cube);
% 
% [X,Y,Z] = mkcube(0.2, 'pose', SE3(0, 0, 1), 'edge');
% cam.mesh(X, Y, Z)
% Tcam = SE3(-1,0,0.5)*SE3.Ry(0.8);
% cam.mesh(X, Y, Z, 'pose', Tcam);
% 
%  theta = [0:500]/100*2*pi;
%  [X,Y,Z] = mkcube(0.2, [], 'edges');
%  for th=theta
%     T_cube = SE3(0, 0, 1.5)*SE3.rpy(th*[1.1 1.2 1.3])
%     cam.mesh( X, Y, Z, 'objpose', T_cube ); drawnow
%  end
%  

%  11.2 Camera Calibration
% P = mkcube(0.2);
% T_unknown = SE3(0.1, 0.2, 1.5) * SE3.rpy(0.1, 0.2, 0.3);
% cam = CentralCamera('focal', 0.015, ...
% 'pixel', 10e-6, 'resolution', [1280 1024], 'noise', 0.05);
% p = cam.project(P, 'objpose', T_unknown);
% C = camcald(P, p) %标定函数 P为笛卡尔空间的坐标，p为图像中的坐标
% 
% null(C)'
% h2e(ans)'
% T_unknown.inv.t'
% est = invcamcal(C)  %分解相机矩阵
% 
% cam.f/cam.rho(2)
% cam.f/cam.rho(2)  %求解行列像素个数
% 
% trprint(T_unknown*est.T)
% hold on; plot_sphere(P, 0.03, 'r')
% trplot(eye(4,4), 'frame', 'T', 'color', 'b', 'length', 0.3)
% est.plot_camera()
% 
% cam = CentralCamera('focal', 0.015, 'pixel', 10e-6, ...
% 'resolution', [1280 1024], 'centre', [640 512]);
% P = mkcube(0.2);
% T_unknown = SE3(0.1, 0.2, 1.5) * SE3.rpy(0.1, 0.2, 0.3);
% p = cam.project(P, 'objpose', T_unknown);
% T_est = cam.estpose(P, p).print
% 
% calib_gui

% 11.3 Wide Field-of-View Imaging
% cam = FishEyeCamera('name', 'fi sheye', ...   %鱼眼镜头
% 'projection', 'equiangular', ...
% 'pixel', 10e-6, ...
% 'resolution', [1280 1024])
% [X,Y,Z] = mkcube(0.2, 'centre', [0.2, 0, 0.3], 'edge');
% cam.mesh(X, Y, Z)
% 
% 
% cam = CatadioptricCamera('name', 'panocam', ...   %射线反射镜
% 'projection', 'equiangular', ...  
% 'maxangle', pi/4, ...
% 'pixel', 10e-6, ...
% 'resolution', [1280 1024])
% [X,Y,Z] = mkcube(1, 'centre', [1, 1, 0.8], 'edge');
% cam.mesh(X, Y, Z)
% 
% cam = SphericalCamera('name', 'spherical') %球镜
% [X,Y,Z] = mkcube(1, 'centre', [2, 3, 1], 'edge');
% cam.mesh(X, Y, Z)

% 11

% L = Plucker([0 0 1], [1 1 1]);
% cam = CentralCamera('default');
% l = cam.project(L)';
% cam.plot(l)
% cam = CentralCamera('default', 'pose', SE3(0.2,0.1, -5)*SE3.Rx(0.2));
% Q = diag([1 1 1 -1]);
% Qs = inv(Q)*det(Q);
% cs = cam.C * Qs * cam.C';
% c = inv(cs)*det(cs);
% det(c(1:2,1:2))

%% Chapter 12 Images and Image Processing
% 12.1 Obtaining an Image
% street = iread('street.png');
% street(200,300)
% 
% %unit8相关运算
% a = uint8(100);
% b = uint8(200);
% a+b;
% 
% %读取与绘制图片
% streetd = iread('street.png', 'double');
% idisp(street)   
% flowers = iread('flowers8.png');
% pix = flowers(276,318,:)
% squeeze(pix)';
% idisp( flowers(:,:,1) );
% seq = iread('seq/*.png');
% 
% %md包含照片其他信息和相机信息
% [im,md]=iread('church.jpg');
% md.DigitalCamera
% 
% %从movie中导入
% cam = Movie('traffic_sequence.mpg');
% cam.size()
% im = cam.grab();
% %播放视频 p365查看相关参数
% while 1
% im = cam.grab;
% if isempty(im) break; end
%  image(im); drawnow
% end
%  
% %从网页获取图片
% cam = AxisWebCamera('http://wc2.dartmouth.edu');
% cam.size();
% im = cam.grab();
% 
% % Images from Maps
% ev = EarthView('key', YOUR_KEY);

% Images from Code
% im = testpattern('rampx', 256, 2);
% im = testpattern('siny', 256, 2);
% im = testpattern('squares', 256, 50, 25);
% im = testpattern('dots', 256, 256, 100);
% canvas = zeros(1000, 1000);
% sq1 = 0.5 * ones(150, 150);
% sq2 = 0.9 * ones(80, 80);
% canvas = ipaste(canvas, sq1, [100 100]);
% canvas = ipaste(canvas, sq2, [300 300]);
% circle = 0.6 * kcircle(120);
% canvas = ipaste(canvas, circle, [600, 200]);
% canvas = iline(canvas, [100 100], [800 800], 0.8);
% 
% %柱状图
% church = iread('church.png', 'grey');
% ihist( church );
% [n,v] = ihist(church);
% [~,x] = peak(n, v);
% [~,x] = peak(n, v,' scale', 25) %比周围X个元素的值都要大的值

% 12.3Monadic Operations
% flowers = iread('flowers8.png');
% church = iread('church.png', 'grey');
% imd = idouble(church);
% im = iint(imd);
% grey = imono(flowers);
% color = icolor(grey); %这是怎么恢复的啊
% color = icolor(grey, [1 0 0]);
% 
% %二值处理
% bright = (church >= 180);
% idisp(bright)
% 
% %直方图均衡化
% im = istretch(church);
% im = inormhist(church);
% ihist(church, 'cdf');
% 
% im = igamm(church, 1/0.45);
% im = igamm(church, 'sRGB');
% idisp( church/64 )

% 12.4Diadic Operations
% subject = iread('greenscreen.jpg', 'double');
% linear = igamm(subject, 'sRGB');
% [r,g] = tristim2cc(linear);
% ihist(g)
% mask = g < 0.45;   %二值分割物体与背景
% idisp(mask) 
% mask3 = icolor( idouble(mask) );
% idisp(mask3 .* subject);
% bg = iread('road.png', 'double'); %背景色
% bg = isamesize(bg, subject);
% idisp(bg .* (1-mask3))  %去除物体显示
% idisp( subject.*mask3 + bg.*(1-mask3) ); %合成显示
% 
% ipixswitch(mask, subject, bg); %图像合成函数，蒙版，物体和背景
% 
% 
% vid = Movie('traffi c_sequence.mpg', 'grey', 'double');

% 12.5 Spatial Operations
%均值处理
% K = ones(21,21) / 21^2;
% mona = iread('monalisa.png', 'double', 'grey');
% idisp( iconvolve(mona, K) );
% % 高斯模糊
% K = kgauss(5);
% idisp( iconvolve(mona, K) );
% idisp( ismooth(mona, 5) )  %模糊函数
% 
% idisp( K );%核函数
% surfl (-15:15, -15:15, K); %help  surfl
% 
% K = kcircle(8, 15);
% 
% %边界检测
% castle = iread('castle.png', 'double', 'grey');
% p = castle(360,:);
% plot(p);
% plot(diff(p))
% 
% K = [0.5 0 -0.5];
% idisp( iconvolve(castle, K), 'invsigned')
% 
% Du = ksobel   %sobel 算子 
% %相关数字图像算法
% idisp( iconvolve(castle, Du), 'invsigned')
% Gu = iconvolve( Du, kgauss(sigma) , 'full');
% edges = icanny(castle, 2);
% 
% mona = iread('monalisa.png', 'double', 'grey');
% T = mona(170:220, 245:295);
% %绝对差之和，均方差之和，相关系数
% sad(T, T)
% ssd(T, T)
% ncc(T, T)
% sad(T, T*0.9)
% ssd(T, T*0.9)
% ncc(T, T*0.9)
% ncc(T, T+0.1)
% zsad(T, T+0.1) %减去均值在进行相关操作
% 
% %相关度检测
% crowd = iread('wheres-wally.png', 'double');
% idisp(crowd)
% wally = iread('wally.png', 'double');
% idisp(wally)
% S = isimilarity(wally, crowd, @zncc);
% idisp(S, 'colormap', 'jet', 'bar')
% 
% idisp(crowd);
% [mx,p] = peak2(S, 1, 'npeaks', 5);
% 
% % 12.6Mathematical Morphology
% eg_morph1
% S = ones(5,5);
% mn = imorph(im, S, 'min');
% morphdemo(im, S, 'min')


% 12.7 Shape Changing
% copy
% mona = iread('monalisa.png');
% [eyes,roi] = iroi(mona);
% idisp(eyes)
% smile = iroi(mona, [265 342; 264 286]);
% %resize
% roof = iread('roof.jpg', 'grey');
% smaller = roof(1:7:end,1:7:end);
% smaller = idecimate(roof, 7);
% bigger = ireplicate( smaller, 7 );
% smoother = ismooth( bigger, 4);
% smaller = iscale(lena, 0.1);%用scale实现缩放
% 
% p = ipyramid( imono(mona) )
% 
% % 缩放
% mona = iread('monalisa.png', 'double', 'grey');
% [Ui,Vi] = imeshgrid(mona);
% [Up,Vp] = imeshgrid(400, 400);
% U = 4*(Up-100); V = 4*(Vp-200);
% little_mona = interp2(Ui, Vi, mona, U, V);
% 
% %旋转pi/6
% R = SO2(pi/6).R; uc = 256; vc = 256;
% U = R(1,1)*(Up-uc) + R(2,1)*(Vp-vc) + uc;
% V = R(1,2)*(Up-uc) + R(2,2)*(Vp-vc) + vc;
% twisted_mona = interp2(Ui, Vi, mona, U, V);
% twisted_mona = irotate(mona, pi/6);
% 
% distorted = iread('Image18.tif', 'double');

%% Chapter 13 Image Feature Extraction
% 13.1 Region Features
% castle = iread('castle.png', 'double');
% idisp(castle >= 0.7)
% ithresh(castle)
% ithresh(castle*0.8)
% ihist(castle);%借助直方图选择阈值
% 
% t = otsu(castle)  %选取阈值的一个方式
% 
% t = niblack(castle, -0.1, 30);
% idisp(castle >= t)   %利用局部阈值进行划分
% 
% [mser,nsets] = imser(castle, 'area', [100 20000]);  %需要编译 %返回值包含照片和类数
% idisp(mser, 'colormap', 'jet') %显示不同的分类结果

% im_targets = iread('yellowtargets.png');
% im_garden = iread('tomato_124.jpg');
% 
% % [cls, cab,resid] = colorkmeans(im_targets, 2, 'ab');   %vl_kmeans函数未能编译成功
% vl_demo_sift_basic
% im = iread('multiblobs.png');
% idisp(im)
% [label, m] = ilabel(im);
% % [label] = ilabel(im);
% idisp(label, 'colormap', jet, 'bar')
% reg3 = (label==3);
% sum(reg3(:))
% [label, m, parents, cls] = ilabel(im);
% parents'  %描述拓扑和等级关系，即该区域被哪些区域包括
% cls'  %不同区域之间的像素值的分类
% targets_label = ilabel(targets_binary);


% im = iread('58060.jpg');
% [label] = igraphseg(im, 1500, 100, 0.5);
% idisp(label, 'colormap', 'jet')

% Bounding Boxes
% sharks = iread('sharks.png');
% [label, m] = ilabel(sharks);
% blob = (label == 2);
% %确定边框
% [v,u] = find(blob);
% umin = min(u)
% umax = max(u)
% vmin = min(v)
% vmax = max(v)
% 
% m00 = mpq(blob, 0, 0) %矩计算
% uc = mpq(blob, 1, 0) / m00
% vc = mpq(blob, 0, 1) / m00
% hold on; plot(uc, vc, 'gx', uc, vc, 'go');
% 
% f = imoments(blob)%直接计算面积，中心，偏角，长短轴比
% fv = iblobs(sharks) %计算各个区域的面积转矩等；功能十分强大的函数
% fv(2). plot_box('g')
% [fv,L] = iblobs(sharks, 'class', 1);
% for i=1:4
%     H(i,:) = humoments(L == fv(i).label);
% end
% 
% fv = iblobs(sharks, 'boundary', 'class', 1);%获取边界
% 
% 
% castle = iread('castle.png');
% words = ocr(castle, [420 300 580 420]);%一个很强大的函数，文字识别
% words.Text
% words.WordConfidences'
% plot_box('matlab', words.WordBoundingBoxes, 'y')

% 13.2Line Features

% im = iread('5points.png', 'double');
% im = testpattern('squares', 256, 256, 128);
% im = irotate(im, -0.3);
% edges = icanny(im);
% h = Hough(edges)
% h.show();
% lines = h.lines()
% %选取距离不同的结果
% h = Hough(edges, 'suppress', 5)
% lines = h.lines()
% 
% idisp(im);
% h.plot('b')
% 
% %复杂图片分割流程
% im = iread('church.png', 'grey', 'double');
% edges = icanny(im);
% h = Hough(edges, 'suppress', 10);
% lines = h.lines();
% idisp(im, 'dark');
% lines(1:10).plot('g');
% lines = lines.seglength(edges);
% k = find( lines.length > 80);
% lines(k).plot('b--')

% 13.3Point Features
%点探测方法
% b1 = iread('building2-1.png', 'grey', 'double');
% C = icorner(b1, 'nfeat', 200);
% idisp(b1, 'dark');
% C.plot('ws');
% Cs = icorner(b1, 'nfeat', 200, 'suppress', 10);
% [C,strength] = icorner(b1, 'nfeat', 200);
% idisp(strength, 'invsigned')
% 
% b2 = iread('building2-2.png', 'grey', 'double');
% C2 = icorner(b2, 'nfeat', 200);   %查询相关帮助文档
% idisp(b2,'dark')
% C2.plot('ws');
% 
% %scale影响
% im = iread('scale-space.png', 'double'); %图片读取有问题
% [G,L,s] = iscalespace(im, 60, 2);
% idisp(L(:,:,5), 'invsigned')
% 
% im = iread('lena.pgm', 'double');%图片读取不了
% % isurf       % Speeded Up Robust Feature%
% sf1 = isurf(b1, 'nfeat', 200)
% idisp(b1, 'dark');
% sf1.plot_scale('g', 'clock')
% hist(sf1.scale, 100);


%% Chapter 14Using Multiple Images
% % 14.1 Feature Correspondence
% im1 = iread('eiffel2-1.jpg', 'mono', 'double');
% im2 = iread('eiffel2-2.jpg', 'mono', 'double');
% 
% %特征点的提取
% hf = icorner(im1, 'nfeat', 200);
% idisp(im1, 'dark'); hf.plot('gs');
% sf = isurf(im1, 'nfeat', 200);
% idisp(im1, 'dark'); sf.plot_scale('g');
% 
% hf(1).descriptor' %0.0805    0.0821    0.0371
% hf(1).distance( hf(2) )%两个特征点之间的欧式距离，但是只有三个特征值
% 
% hf = icorner(im1, 'nfeat', 200, 'color', 'patch', 5)
% hf(1).ncc( hf(2) )  %描述点之间是否匹配的算法，值为[-1,1]，>0.8视为可信
% %运算成本较高
% 
% %寻求匹配点算法
% sf1 = isurf(im1);
% sf2 = isurf(im2)
% m = sf1.match(sf2)  
% m(1:5)
% 
% idisp({im1, im2}, 'dark')
% m.subset(100).plot('w') %连接对应点
% %求解对应点的编号
% [m,corresp] = sf1.match(sf2);
% corresp(:,1:5)  
% 
% %绘制直方图，确定大概有多少比较靠近
% m2 = sf1.match(sf2, 'all');
% histogram(m2.distance, 'Normalization', 'cdf')
% 
% mm = sf1.match(sf2, 'thresh', 0.05); %距离小于0.05距离值的匹配点
% %选取匹配度最高的N个点
% N=20
% mm = sf1.match(sf2, 'top', N);   %问题：如何实现匹配算法的
% 
% 
% % 14.2Geometry of Multiple Views
% %确定相机1和相机2的相关参数
% T1 = SE3(-0.1, 0, 0) * SE3.Ry(0.4);
% cam1 = CentralCamera('name', 'camera 1', 'default', ...	
% 'focal', 0.002, 'pose', T1);
% T2 = SE3(0.1, 0,0)*SE3.Ry(-0.4);
% cam2 = CentralCamera('name', 'camera 2', 'default', ...	
% 'focal', 0.002, 'pose', T2);
% %绘制相机位置
% axis([-0.5 0.5 -0.5 0.5 0 1]);
% cam1.plot_camera('color', 'b', 'label')
% cam2.plot_camera('color', 'r', 'label')
% 
% P=[0.5 0.1 0.8]';
% plot_sphere(P, 0.03, 'b');
% %点在相机投影
% p1 = cam1.plot(P)
% p2 = cam2.plot(P)
% %相机互相投影
% cam1.hold
% e1 = cam1.plot( cam2.centre, 'Marker', 'd',	...
% 'MarkerFaceColor', 'k');
% cam2.hold
% e2 = cam2.plot( cam1.centre, 'Marker', 'd',	...
% 'MarkerFaceColor', 'k');
% 
% F = cam1.F( cam2 ) %为什么这么求F
% e2h(p2)' * F * e2h(p1) %欧式坐标与射影坐标坐转换
% rank(F)
% null(F)'
% e1 = h2e(ans)'
% null(F')
% e2 = h2e(ans)'
% 
% cam2.plot_epiline(F, p1, 'r')
% cam1.plot_epiline(F', p2, 'r');
% 
% E = cam1.E(F)
% sol = cam1.invE(E)
% inv(cam1.T) * cam2.T %相机1相对于相机2的位置
% 
% Q = [0 0 10]'; %利用已知点确定相关参数
% cam1.project(Q)'
% cam1.move(sol(1).T).project(Q)'
% cam1.move(sol(2).T).project(Q)'
% sol = cam1.invE(E, Q)
% 
% 
% P = SE3(-1, -1, 2)*(2 *rand(3,20) );
% p1 = cam1.project(P);
% p2 = cam2.project(P);
% F = fmatrix(p1, p2)            %这个函数的功能是啥？
% cam2.plot(P);
% 
% p2(:,[8 7]) = p2(:,[7 8]); %错误点示范
% epidist(F, p1(:,1), p2(:,1))
% epidist(F, p1(:,7), p2(:,7))
% 
% [F,in,r] = ransac(@fmatrix, [p1; p2], 1e-6, 'verbose'); %in代表正确的点对应的位置
% F = m.ransac(@fmatrix, 1e-4, 'verbose') 
% m.show
% m(1:5) %对应点，距离
% %两类点绘制
% idisp({im1, im2});
% m.inlier.subset(100).plot('g')
% idisp({im1, im2});
% m.outlier.subset(100).plot('r')
% 
% cam = CentralCamera('image', im1)
% cam.plot_epiline(F', m.inlier.subset(20).p2, 'g');
% 
% h2e( null(F))
% cam.plot(ans, 'bo')
% 
% Tgrid = SE3(0,0,1)*SE3.Rx(0.1)*SE3.Ry(0.2);
%c
% p1 = cam1.plot(P, 'o');
% p2 = cam2.plot(P, 'o');
% H = homography(p1, p2)
% p2b = homtrans(H, p1);
% cam2.hold()
% cam2.plot(p2b, '+')
% p1b = homtrans(inv(H), p2);
% Q = [
% -0.2302 -0.0545 0.2537
% 0.3287 0.4523 0.6024
% 0.4000 0.5000 0.6000 ];
% axis([-1 1 -1 1 0 2])
% plot_sphere(P, 0.05, 'b')
% plot_sphere(Q, 0.05, 'r')
% cam1.plot_camera('color', 'b', 'label')
% cam2.plot_camera('color', 'r', 'label')
% p1 = cam1.plot([P Q], 'o');
% p2 = cam2.plot([P Q], 'o');
% p2h = homtrans(H, p1);
% cam2.plot(p2h, '+')
% 
% colnorm( homtrans(H, p1)-p2 ) %判断是否在一个平面。
% [H,in] = ransac(@homography, [p1; p2], 0.1)
% cam1.invH(H)
% inv(T1)*T2
% inv(T1)* Tgrid
% 
% 
% im1 = iread('walls-l.jpg', 'double', 'reduce', 2);
% im2 = iread('walls-r.jpg', 'double', 'reduce', 2);
% sf1 = isurf(im1);
% sf2 = isurf(im2);
% m = sf1.match(sf2, 'top', 1000)
% [H,r] = m.ransac(@homography, 4)
% m.show
% idisp(im1)
% plot_point(m.inlier.p1, 'ys') %绘制匹配点
% m = m.outlier

% % 14.3 Stereo Vision
% [F,r] = m.ransac(@fmatrix,1e-4, 'verbose');
% cam = CentralCamera('image', im1);
% cam.plot_epiline(F', m.inlier.subset(40).p2, 'y');
% [~,md] = iread('walls-l.jpg'); %md is a structure of text strings that contains various characteristics of the image its metadat
% f = md.DigitalCamera.FocalLength
% md.Model
% cam = CentralCamera('image', im1, 'focal', f/1000, ...	
% 'pixel', 2*1.5e-6)
% E = cam.E(F)
% T = cam.invE(E, [0,0,10]')
% T.torpy('yxz', 'deg')
% t = T.t;
% %绘制两条射线与交点
% r1 = cam.ray(m(1).p1)
% r2 = cam.move(T).ray(m(1).p2)
% [P,e] = r1.intersect(r2);
% 
% m2 = m.inlier.subset(100);
% r1 = cam.ray( m2.p1 );
% r2 = cam.move(T).ray( m2.p2 );
% [P,e] = r1.intersect(r2);
% z = P(3,:)
% idisp(im1)
% plot_point(m2.p1, 'y+', 'textcolor', 'y', 'printf', {'%.1f', z});

% 14.3.2 Dense Stereo Matching
% L = iread('rocks2-l.png', 'reduce', 2);
% R = iread('rocks2-r.png', 'reduce', 2);
% d = istereo(L, R, [40, 90], 3);
% idisp(d, 'bar')
% [d,sim,DSI] = istereo(L, R, [40 90], 3);%p484测量水平移动像素的相似性
% plot( squeeze(DSI(439,138,:)), 'o-');
% idisp(sim)
% ipixswitch(sim<0.7, 'yellow', d/90);
% ihist(sim(isfinite(sim)), 'normcdf');
% slice(DSI, [], [100 200 300 400 500], [])
% shading interp; colorbar
% 
% [di,sim,peak] = istereo(L, R, [40 90], 3, 'interp');
% peak
% % 14.3.4Cleaning up and Reconstruction
% status = ones(size(d));
% [U,V] = imeshgrid(L);
% status(isnan(d)) = 5; % search template off the edge
% status(U<=90) = 2; % no overlap
% status(sim<0.8) = 3; % weak match
% status(peak.A>=-0.1) = 4; % broad peak
% idisp(status)
% colormap( colorname({'lightgreen', 'cyan', 'blue', 'orange', 'red'}) )
% sum(status(:) == 1) / prod(size(status)) * 100
% di(status>1) = NaN;
% ipixswitch(isnan(di), 'red', di/90);
% 
% di = di + 274; %为什么是274
% [U,V] = imeshgrid(L);
% u0 = size(L,2)/2; v0 = size(L,1)/2;
% b = 0.160;
% X = b*(U-u0) ./ di; Y = b*(V-v0) ./ di; Z = 3740 * b ./ di;
% surf(Z)
% shading interp; view(-150, 75)
% set(gca,'ZDir', 'reverse'); set(gca,'XDir', 'reverse')
% colormap(flipud(hot))
% 
% dimf = irank(di, 41, ones(9,9));%filter设置
% di = ipixswitch(isnan(di), dimf, di);
% X = b*(U-u0) ./ di; Y = b*(V-v0) ./ di; Z = 3740 * b ./ di;
% 
% Lcolor = iread('rocks2-l.png');
% surface(X, Y, Z, Lcolor, 'FaceColor', 'texturemap', ...
% 'EdgeColor', 'none', 'CDataMapping', 'direct')
% xyzlabel
% set(gca,'ZDir', 'reverse'); set(gca,'XDir', 'reverse')
% 
% anaglyph(L, R, 'rc') %立体视觉
% L = iread('walls-l.jpg', 'mono', 'double', 'reduce', 2);
% R = iread('walls-r.jpg', 'mono', 'double', 'reduce', 2);
% sL = isurf(L);
% sR = isurf(R);
% m = sL.match(sR, 'top', 1000);
% F = m.ransac(@fmatrix,1e-4, 'verbose');
% [Lr,Rr] = irectify(F, m, L, R);
% stdisp(Lr, Rr)
% d = istereo(Lr, Rr, [180 530], 7, 'interp');

% 14.4 Bundle Adjustment
% im1 = iread('walls-l.jpg', 'double', 'reduce', 2);
% cam = CentralCamera('image', im1);
%  [~,md] = iread('walls-l.jpg');
% f = md.DigitalCamera.FocalLength;
% cam = CentralCamera('image', im1, 'focal', f/1000, ...	
% 'pixel', 2*1.5e-6)
% P = SE3(-1, -1, 2)*(2 *rand(3,20) );
% p1 = cam.project(P);
% p2 = cam.move(T).project(P);
% e = colnorm( [p1-m2.p1 p2-m2.p2] );

% 14.5 Point Clouds
% T = SE3(1,2,3) * SE3.rpy(0.3, 0.4, 0.5);
% P = mkgrid(10, 1, 'pose', T);
% P = P + 0.02*randn(size(P));
% x0 = mean(P')
% P = bsxfun(@minus, P, x0');%减去均值
% J = P*P'
% [x,lambda] = eig(J); %特征值与特征向量
% diag(lambda)'
% n = x(:,1)' %第一组特征向量,为啥是平面法向量
% T. SO3.a'   %z轴方向，为啥会成为这个平面的法向量
% 
% % P505
% load bunny
% M = bunny;
% T_unknown = SE3(0.2, 0.2, 0.1) * SE3.rpy(0.2, 0.3, 0.4);
% D = T_unknown * M;
% corresp = closest(D, M);
% [T,d] = icp(M, D, 'plot'); %如果需要查看相关函数，匹配
% trprint(T, 'rpy', 'radian')
% %验证函数的鲁棒性
% D(:,randi(numcols(D), 40,1)) = [];
% D = [D 0.1*rand(3,20)+0.1];
% D = D + 0.01*randn(size(D));
% [T,d] = icp(M, D, 'plot', 'distthresh', 3);
% trprint(T, 'rpy', 'radian')

%14.6 Structured Light

% 14.7 Applications
% im = iread('notre-dame.jpg', 'double');
% idisp(im)
% p1 = ginput(4)'
% plot_poly(p1, 'wo', 'fill', 'b', 'alpha', 0.2);
% %绘制矩形框
% mn = min(p1');
% mx = max(p1');
% p2 = [mn(1) mx(2); mn(1) mn(2); mx(1) mn(2); mx(1) mx(2)]';
% plot_poly(p2, 'k', 'fill', 'r', 'alpha', 0.2)
% H = homography(p1, p2) %转换矩阵，随机框转换为矩形框
% homwarp(H, im, 'full')
% 
% [~,md] = iread('notre-dame.jpg', 'double');
% f = md.DigitalCamera.FocalLength  %单位毫米
% cam = CentralCamera('image', im, 'focal', f/1000, ...
% 'sensor', [7.18e-3,5.32e-3])
% sol = cam.invH(H, 'verbose');
% tr2rpy(sol(2).T, 'deg', 'camera') %求取变换矩阵的RPY
% 
% % Mosaicing 拼接照片
% im1 = iread('mosaic/aerial2-1.png', 'double', 'grey');
% im2 = iread('mosaic/aerial2-2.png', 'double', 'grey');
% composite = zeros(2000,2000);
% composite = ipaste(composite, im1, [1 1]);
% %特征点匹配
% f1 = isurf(im1)
% f2 = isurf(im2)
% m = f1.match(f2);
% [H,in] = m.ransac(@homography, 0.2)
% [tile,t] = homwarp(inv(H), im2, 'full', 'extrapval', 0); %调整函数
% composite = ipaste(composite, tile, t, 'add');%拼接函数
% [tile,t] = homwarp(inv(H), im2, 'full', 'extrapval', NaN);
% composite = ipaste(composite, tile, t, 'mean');
% % Image Matching and Retrieval
% images = iread('campus/*.jpg', 'mono');
% sf = isurf(images, 'thresh', 0);  %各个照片的特征点
% bag = BagOfWords(sf, 2000) %标记
% w = bag.words(259);
% [word,f] = bag.wordfreq() ; %word is a vector containing all unique words and f are their corresponding
% %frequencies.
% bar( sort(f, 'descend') )
% bag.remove_stop(50)
% M = bag.wordvector;%词向量
% S = bag.similarity(bag) %利用词向量之间的夹角值可以判断其相似性
% % Visual Odometry vodemo
% left = iread('bridge-l/*.png', 'roi', [20 750; 20 440]);
% ianimate(left, 'fps', 10);%播放照片
% c = icorner(left, 'nfeat', 200, 'patch', 7);
% ianimate(left, c, 'fps', 10);
% right = iread('bridge-r/*.png', 'roi', [20 750; 20 440]);

%这一章内容再看一遍2019.8.14
%% Part V Robotics, Vision and Control
% % Chapter15 Vision-Based Control
% 15.1 Position-Based Visual Servoing
% cam = CentralCamera('default');%相机
% P = mkgrid( 2, 0.5, 'pose', SE3(0,0,3) ); %目标点
% p = cam.plot(P, 'pose', T_C); %T_C为设置
% C_Te_G = cam. estpose(P, p);%相机与目标相对位姿
% T_delta = C_Te_G * inv(Cd_T_G);
% T_C = T_C .* T_delta.interp(lambda); %lambda未设置？ lambda=0.5

%15.2
% cam = CentralCamera('default');
% P = [1 1 5]';
% p0 = cam.project( P )
% px = cam.project( P, 'pose', SE3(0.1,0,0) )
% %相机不同位置的偏移量
% ( px - p0 ) / 0.1
% ( cam.project( P, 'pose', SE3(0, 0, 0.1) ) - p0 ) / 0.1
% ( cam.project( P, 'pose', SE3.Rx(0.1) ) - p0 ) / 0.1
% 
% %图像雅克比矩阵，2X6
% J = cam.visjac_p([672; 672], 5)
% %相机不同速度是像素的的移动
% cam.flowfield( [1 0 0 0 0 0] );
% cam.flowfield( [0 0 1 0 0 0] );
% cam.flowfield( [0 0 0 0 0 1] );
% %投影到中心点，距离1m的图像雅克比
% cam.visjac_p(cam.pp', 1)
% 
% %不同的焦距下旋转的影响
% cam.f = 20e-3;
% cam.flowfield( [0 0 0 0 1 0] );
% cam.f = 4e-3;
% cam.flowfield( [0 0 0 0 1 0] );
% %求取零空间
% J = cam.visjac_p(cam.pp', 1);
% null(J)
% 
% cam = CentralCamera('default');
% P = mkgrid( 2, 0.5, 'pose', SE3(0,0,3) );
% pd = bsxfun(@plus, 200*[-1 -1 1 1; -1 1 1 -1], cam.pp');
% T_C=SE3(0,0.3,0.4);
% p = cam.plot(P, 'pose', T_C);
% e = pd - p;
% J = cam.visjac_p( p, 1 );
% lambda=0.1;
% v = lambda * pinv(J) * e(:);
% T_C = T_C .* delta2tr(v);
% 
% T_C0 = SE3(1,1,-3)*SE3.Rz(0.6);%初始姿态
% ibvs = IBVS(cam, 'pose0', T_C0, 'pstar', pd);
% ibvs.run();
% ibvs.plot_p();
% ibvs.plot_vel();
% ibvs.plot_camera();
% 
% sl_ibvs         %Simulink仿真
% 
% %考虑众多参数不确定时？
% ibvs = IBVS(cam, 'pose0', T_C0, 'pstar', pd,'depth', 1)
% ibvs.run(50)
% ibvs = IBVS(cam, 'pose0', T_C0, 'pstar', pd,'depth', 10)
% ibvs.run(50)

% 15.3 Using Other Image Features
% 15.3.1Line Features
% P = circle([0 0 3], 0.5, 'n', 3); %获取点
% ibvs = IBVS_l(cam, 'example');
% ibvs.run()
% 
% P = circle([0 0 3], 0.5, 'n', 10);
% p = cam.project(P, 'pose', TC);
% pn = cam.normalized( p);
% x = pn(1,:); y = pn(2,:);
% a = [y.^2; -2*x.*y; 2*x; 2*y; ones(1,numcols(x))]';
% b = -(x.^2)';
% E = a\b;

% Chapter  16 Advanced Visual Servoing
% 16.1 XY/Z-Partitioned IBVS
% sl_partitioned %可以看看
% slcamera(cam, u)
% % 16.2 IBVS Using Polar Coordinates
% % 16.3 IBVS for a Spherical Camera
% % 16.4 Applications
% sl_arm_ibvs %多看看
% sl_omni_vs %p574
% sl_drivepose_vs
% sl_quadrotor_vs


%% 移动机器人
% % Part II Mobile Robots
%% Chapter4 Mobile Robot Vehicles
% 4.1 Wheeled Mobile Robots
% 4.1.1 Car-Like Mobile Robots
% sl_lanechange
% t = out.get('t'); q = out.get('y');
% mplot(t, q)
% %点
% sl_drivepoint
% xg = [5 5]; %设定目标值
% x0 = [8 5 pi/2];%设定初始值
% r = sim('sl_drivepoint');
% %线追踪
% sl_driveline
% L = [1 -2 4];
% x0 = [8 5 pi/2];
% r = sim( 'sl_driveline');
% %轨迹
% sl_pursuit
% r = sim('sl_pursuit')
% %位姿
% xg = [5 5 pi/2];
% x0 = [9 5 0];
% sl_drivepose
% % 4.2 Flying Robots
% sl_quadrotor %四旋翼
% mdl_quadrotor

% 4.3 Advanced Topics


%% Chapter 5 Navigation
% % 5.1 Reactive Navigation
% % sl_braitenberg
% %  function sensor = sensorfield(x, y)
% %  xc = 60; yc = 90;
% %  sensor = 200./((x-xc).^2 + (y-yc).^2 + 200);
% %  end
%  
%  load house
% %  place
%  bug = Bug2(house);
%  bug.plot();
%  bug.query(place.br3, place.kitchen, 'animate');
% %  p = bug.query(place.br3, place.kitchen);
%  p = bug.query([], place.kitchen,'animate') ;
%  
% %  5.2 Map-Based Planning
% dx = DXform(house); %距离，得到一个类
% dx.plan(place.kitchen) %需要安装视觉工具包
% dx.query(place.br3, 'animate');
% p = dx.query(place.br3);
% dx.plot3d(p) %3d显示距离
% 
% ds = Dstar(house); %D*方法
% c = ds.costmap();
% ds.plan(place.kitchen);
% ds.niter
% ds.query(place.br3);
% ds.modify_cost( [300,325; 115,125], 5 ); %更改cost函数
% ds.niter
% ds.query(place.br3);
% 
% % Roadmap Methods
% free = 1 - house;
% free(1,:) = 0; free(end,:) = 0;
% free(:,1) = 0; free(:,end) = 0;
% skeleton = ithin(free); %thinning收缩，图片处理
% 
% % Probabilistic Roadmap Method
% %相比于上面的方法，该方法只用绘制一次图，但是得到的路线并不一定最近
% prm = PRM(house)
% prm.plan('npoints', 150)
% prm
% prm.plot()
% prm.query(place.br3, place.kitchen)
% prm.plot()
% p = prm.query(place.br3, place.kitchen);
% 
% % Lattice Planner
% %考虑到车辆无法侧移
% lp = Lattice();
% lp.plan('iterations', 2)
% lp.plot()
% lp.plan('iterations', 8)
% lp.plot()
% lp.query( [1 2 pi/2], [2 -2 0] );
% lp.plot
% 
% lp.plan('cost', [1 10 10])%设置三个路径方法（直走，左右圆弧）的cost
% lp.query(start, goal);
% lp.plot()
% 
% % Rapidly-Exploring Random Tree (RRT)
% %这种搜索是否会陷入局部最优，甚至局部最优和全局最优冲突，全局随机树可以避免？
% car = Bicycle('steermax', 0.5);
% rrt = RRT(car, 'npoints', 1000)
% rrt.plan();
% rrt.plot();
%  load road
% rrt = RRT(car, road, 'npoints', 1000, 'root', [50 22 0], 'simtime', 4) %没有road数据
% rrt.plan();
% p = rrt.query([40 45 0], [50 22 0]);

%% Chapter 6
% 6.1 Dead Reckoning
% V = diag([0.02, 0.5*pi/180].^2);
% veh= Bicycle('covar', V)            %研究自行车很重要的一个类
% odo = veh.step(1, 0.3)
% veh.f([0 0 0], odo)%实际测量
% veh.add_driver( RandomPath(10) )%随机生成点
% veh.run()
% 
% veh.Fx( [0,0,0], [0.5, 0.1] )
% P0 = diag([0.005, 0.005, 0.001].^2); %起始点
% ekf = EKF(veh, V, P0);    %扩展科尔曼滤波
% ekf.run(1000);
% veh.plot_xy()
%  hold on
% ekf.plot_xy('r') %查看估计偏差
% ekf.plot_ellipse('g')

% 6.2 Localizing with a Map  %没大看懂
% map = LandmarkMap(20, 10)
% W = diag([0.1, 1*pi/180].^2);
% sensor = RangeBearingSensor(veh, map, 'covar', W)
% [z,i] = sensor.reading()
% landmark(17)

%相比于之前的估计
% map = LandmarkMap(20);
%  veh = Bicycle('covar', V);
% veh.add_driver( RandomPath(map.dim) );
% sensor = RangeBearingSensor(veh, map, 'covar', W, 'angle',...
% [-pi/2 pi/2], 'range', 4, 'animate');
% ekf = EKF(veh, V, P0, sensor, W, map);
% ekf.run(1000)
% map.plot()
% veh.plot_xy();
% ekf.plot_xy('r');
% ekf.plot_ellipse('k')
% 
% % 6.3 Creating a Map
% map = LandmarkMap(20);
% veh = Bicycle(); % error free vehicle
% veh.add_driver( RandomPath(map.dim) );
% W = diag([0.1, 1*pi/180].^2);
% sensor = RangeBearingSensor(veh, map, 'covar', W);
% ekf = EKF(veh, [], [], sensor, W, []);
% ekf.run(1000);

% 6.4 Localization and Mapping
%  P0 = diag([.01, .01, 0.005].^2);
%  V = diag([0.02, 0.5*pi/180].^2);
%   W = diag([0.1, 1*pi/180].^2);
% map = LandmarkMap(20);
% veh = Bicycle('covar', V);
% veh.add_driver( RandomPath(map.dim) );
% sensor = RangeBearingSensor(veh, map, 'covar', W);
% ekf = EKF(veh, V, P0, sensor, W, []);
% 
% ekf.run(1000)
% 
% map.plot();
% ekf.plot_map('g');
% ekf.plot_xy('r');
% veh.plot_xy('b');

% 6.5 Rao-Blackwellized SLAM
% 6.6 Pose Graph SLAM
%求解雅克比矩阵
%需要用两个位姿才能估计出待测点位置，不过好像一个也可以
% syms xi yi ti xj yj tj xm ym tm assume real
% xi_e = inv( SE2(xm, ym, tm) ) * inv( SE2(xi, yi, ti) ) * SE2(xj, yj, tj);
% fk = simplify(xi_e.xyt);
% jacobian ( fk, [xi yi ti] );
% Ai = simplify (ans)
% pg = PoseGraph('pg1.g2o')
% pg.optimize('animate')%优化位姿估计结果，但是是如何做到的
% %another example
% pg = PoseGraph('killian-small.toro');
% pg.plot()
% 
% % 6.7 Sequential Monte-Carlo Localization
% map = LandmarkMap(20);
% V= diag([0.1, 1*pi/180].^2);
% veh = Bicycle('covar', V);
% veh.add_driver( RandomPath(10) );
% W= diag([0.005, 0.5*pi/180].^2);
% sensor = RangeBearingSensor(veh, map, 'covar', W);
% 
% Q = diag([0.1, 0.1, 1*pi/180]).^2;
% L = diag([0.1 0.1]);
% pf = ParticleFilter(veh, sensor, Q, L, 1000);
% pf.run(1000);
% 
% map.plot();
% veh.plot_xy('b');
% pf.plot_xy('r');
% plot(pf.std(1:100,:))
% pf.plot_pdf()
% 
% % 6.8 Application: Scanning Laser Rangefinder
% pg = PoseGraph('killian.g2o', 'laser');
% [r, theta] = pg.scan(2580);
% polar(theta, r)
% [x,y] = pol2cart (theta, r); %极坐标转换为笛卡尔
% plot (x, y, '.')
% p2580 = pg.scanxy(2580);
% p2581 = pg.scanxy(2581);
% 
% T = icp( p2581, p2580, 'verbose' , 'T0', transl2(0.5, 0), 'distthresh', 3)%函数必须在mex目录下运行
% pg.time(2581)-pg.time(2580)%两个测量之间的时间
% 
% %map
% pg.scanmap()
% pg.plot_occgrid()