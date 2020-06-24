%基于POE模型下牛顿欧拉迭代算法的机器人逆动力学求解
%目前只支持数字解
%Tor=idyna_poe_dh(robot,q,qd,qdd,gravity,f_end)
%输入robot为机器人模型，包含有运动学参数和动力学参数
%输入q为关节角度（1XN），qd为关节角加速度,qdd为对应关节角加速度
%gravity为重力加速度，f_end为末端执行器施加载荷，
%这两个在不输入时会选用默认值
%输出Torque为关节扭矩
%参考文献为：熊有伦等.机器人学建模规划与控制,7.8节:2018
%A Lie Group Formulation of Robot Dynamics

function Tor=idyna_poe_mdh(robot,q,qd,qdd,gravity,f_end)
%% 设置和获取基本参数
%确保输入角度，角速度，角加速度为行向量
q=q(:)';
qd=qd(:)';
qdd=qdd(:)';
n=robot.n; %机器人连杆数目
z0=[0,0,1]'; %基坐标系，z0的方向
grav= robot.gravity; %重力加速度，可以设置
fend = zeros(6, 1);%末端施加的载荷,W=[fx Fy Fz Mx My Mz]';
debug1=1;%设置调试参数
debug2=0;
if nargin>4%判断是否输入重力加速度
    grav=gravity;
end
if nargin==6%判断外部是否施加载荷
    fend=f_end;
end
%% 判断输入是否符合要求
if numcols(q)~=n||numcols(qd)~=n||numcols(qd)~=n||numrows(q)~=1||...
        numrows(qd)~=1||numrows(qdd)~= 1||length(fend)~=6
    error('输入数据有误')
end
 %% 求解 
 offset=robot.offset;
 n=robot.n;
 V0=zeros(6,1);%设置初始速度
 Vd0=zeros(6,1);%设置初始加速度（此书参考文献中均存在疑问，因为加速度没有被考虑其中）
 
 %% 求解T0和惯性矩阵
 T0={};
for i=1:n
    T0{i}=robot.links(i).A(offset(i)).T;
    M{i}=m_matrix(robot.links(i).m,robot.links(i).r,robot.links(i).I);
    if strcmp(robot.links(i).type,'R')
        twist(i,:)=[0 0 0 0 0 1];
    else
        twist(i,:)=[0 0 1 0 0 0];
    end
end

%% 正向递推
for i=1:n
    T{i}=T0{i}*exp_se(twist(i,:),q(i));
    if i==1
        V{i}=Adjoint_v(inv(T{i}))*V0+twist(i,:)'*qd(i);
        Vd{i}=twist(i,:)'*qdd(i)+Adjoint_v(inv(T{i}))*Vd0-ad_v(twist(i,:)'*qd(i))*(Adjoint_v(inv(T{i}))*V0);
    else
        V{i}=Adjoint_v(inv(T{i}))*V{i-1}+twist(i,:)'*qd(i);
        Vd{i}=twist(i,:)'*qdd(i)+Adjoint_v(inv(T{i}))*Vd{i-1}-ad_v(twist(i,:)'*qd(i))*(Adjoint_v(inv(T{i}))*V{i-1});
    end
end

%% 逆向递推
for i=n:-1:1
    if i==n
        F{i}=Adjoint_v(inv(eye(4)))'*Fend+M{i}*Vd{i}-ad_v(V{i})'*M{i}*V{i};
    else
        F{i}=Adjoint_v(inv(T{i+1}))'*F{i+1}+M{i}*Vd{i}-ad_v(V{i})'*M{i}*V{i};
    end
    tor(i)=twist(i,:)*F{i};
end

%% 惯性矩阵函数
    function  M=m_matrix(m,pci,I)
        M(1:3,1:3)=m*eye(3);
        M(1:3,4:6)=-m*skew_poe(pci);
        M(4:6,1:3)=m*skew_poe(pci);
        M(4:6,4:6)=I-m*skew_poe(pci)*skew_poe(pci);
    end
end

