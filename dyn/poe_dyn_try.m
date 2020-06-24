%ָ�����㷨�������˶���ѧ

clc
clear
mdl_puma560
% p560A
%% ָ����ţ��ŷ��
% % p560_poe=dh2poe(p560); %��DH����ģ��ת��ΪPOEģ��
% % twist=p560_poe.twist; %��ȡ�����ؽڵ�����
% offset = p560.offset; %��ȡ���˳�ʼƫ�ƽǶ�
% n=p560.n;             %��ȡ���˸���
% robot=p560;
% %�����ʼ����
% V0=zeros(6,1);
% Vd0=[0 0 9.8100  0 0  0]';
% Fend=zeros(6,1);
% 
% q=qn;
% qd=qr;
% qdd=qs;
% 
% %��Ⲣ�洢T0�;ֲ�twist
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
% %�����������ٶȺͼ��ٶ�
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

%% ָ��������������
robot=p560;
robot_poe=dh2poe(robot);
n=robot.n;
offset=robot_poe.offset;
twist=robot_poe.twist;

syms q;
%���������ز���
%�˴�����������ĺ͹����������ǰһ��
J={};
for i=1:n
    T0{i}=robot.links(i).A(offset(i)).T;
    if i==1
        T{i}=T0{i};
    else
        T{i}=T{i-1}*T0{i};
    end
    Tc{i}=T{i};
    Tc{i}(1:3,4)=Tc{i}(1:3,4)+robot.links(i).r';%����ʼ���ĵ�λ�˾���
    T_temp=Tc{i};
    J{i}=zeros(6,6);
    for j=1:i%����ط���Ҫ����һ��
        if j~=1
            T_temp=inv(T{i-1})*T_temp;
        end
        Lj=Adjoint_v(T_temp)*twist(i,:)';
        J{i}(:,j)=Lj;
    end
    M{i}=m_matrix(robot.links(i).m,robot.links(i).I);
end

%% ���Ծ�����
function  M=m_matrix(m,I)
M=zeros(6,6);
M(1:3,1:3)=m*eye(3);
M(4:6,4:6)=I;
end



