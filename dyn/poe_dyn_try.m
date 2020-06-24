%指数积算法求解机器人动力学

clc
clear
mdl_puma560
% p560A
%% 指数积牛顿欧拉
% % p560_poe=dh2poe(p560); %将DH参数模型转换为POE模型
% % twist=p560_poe.twist; %获取各个关节的旋量
% offset = p560.offset; %获取连杆初始偏移角度
% n=p560.n;             %获取连杆个数
% robot=p560;
% %定义初始变量
% V0=zeros(6,1);
% Vd0=[0 0 9.8100  0 0  0]';
% Fend=zeros(6,1);
% 
% q=qn;
% qd=qr;
% qdd=qs;
% 
% %求解并存储T0和局部twist
% T0={};
% for i=1:n
%     T0{i}=robot.links(i).A(offset(i)).T;
%     if strcmp(robot.links(i).type,'R')
%         twist(i,:)=[0 0 0 0 0 1];
%     else
%         twist(i,:)=[0 0 1 0 0 0];
%     end
% end
% 
% %正向递推求解速度和加速度
% for i=1:n
% %     T{i}=T0{i}*exp_se(twist(i,:),q(i));
%     T{i}=exp_se(twist(i,:),q(i))*T0{i};
%     M{i}=m_matrix(robot.links(i).m,robot.links(i).r,robot.links(i).I);
% %     if i==1
% %         V{i}=Adjoint_v(inv(T{i}))*V0+twist(i,:)'*qd(i);
% %         Vd{i}=twist(i,:)'*qdd(i)+Adjoint_v(inv(T{i}))*Vd0-ad_v(twist(i,:)'*qd(i))*(Adjoint_v(inv(T{i}))*V0);
% %     else
% %         V{i}=Adjoint_v(inv(T{i}))*V{i-1}+twist(i,:)'*qd(i);
% %         Vd{i}=twist(i,:)'*qdd(i)+Adjoint_v(inv(T{i}))*Vd{i-1}-ad_v(twist(i,:)'*qd(i))*(Adjoint_v(inv(T{i}))*V{i-1});
% %     end
%     if i==1
%         V{i}=Adjoint_v(eye(4))*V0+twist(i,:)'*qd(i);
%         Vd{i}=twist(i,:)'*qdd(i)+Adjoint_v(inv(T{i}))*Vd0-ad_v(twist(i,:)'*qd(i))*(Adjoint_v(inv(T{i}))*V0);
%     else
%         V{i}=Adjoint_v(inv(T{i-1}))*V{i-1}+twist(i,:)'*qd(i);
%         Vd{i}=twist(i,:)'*qdd(i)+Adjoint_v(inv(T{i}))*Vd{i-1}-ad_v(twist(i,:)'*qd(i))*(Adjoint_v(inv(T{i}))*V{i-1});
%     end
% end
% 
% tor=zeros(6,1);
% for i=n:-1:1
%     if i==n
%         F{i}=Adjoint_v(inv(eye(4)))'*Fend+M{i}*Vd{i}-ad_v(V{i})'*M{i}*V{i};
%     else
%         F{i}=Adjoint_v(inv(T{i}))'*F{i+1}+M{i}*Vd{i}-ad_v(V{i})'*M{i}*V{i};
%     end
%     tor(i)=twist(i,:)*F{i};
% end
% 
% 
% idyna_nr_dh(p560,q,qd,qdd)
% tor
% 
% function  M=m_matrix(m,pci,I)
% M(1:3,1:3)=m*eye(3);
% M(1:3,4:6)=-m*skew_poe(pci);
% M(4:6,1:3)=m*skew_poe(pci);
% M(4:6,4:6)=I-m*skew_poe(pci)*skew_poe(pci);
% end

%% 指数积的拉格朗日
robot=p560;
robot_poe=dh2poe(robot);
n=robot.n;
offset=robot_poe.offset;
twist=robot_poe.twist;

syms q;
%迭代求解相关参数
%此处假设的是质心和惯量是相对于前一个
J={};
for i=1:n
    T0{i}=robot.links(i).A(offset(i)).T;
    if i==1
        T{i}=T0{i};
    else
        T{i}=T{i-1}*T0{i};
    end
    Tc{i}=T{i};
    Tc{i}(1:3,4)=Tc{i}(1:3,4)+robot.links(i).r';%求解初始质心的位姿矩阵
    T_temp=Tc{i};
    J{i}=zeros(6,6);
    for j=1:i%这个地方需要考虑一下
        if j~=1
            T_temp=inv(T{i-1})*T_temp;
        end
        Lj=Adjoint_v(T_temp)*twist(i,:)';
        J{i}(:,j)=Lj;
    end
    M{i}=m_matrix(robot.links(i).m,robot.links(i).I);
end

%% 惯性矩阵函数
function  M=m_matrix(m,I)
M=zeros(6,6);
M(1:3,1:3)=m*eye(3);
M(4:6,4:6)=I;
end



