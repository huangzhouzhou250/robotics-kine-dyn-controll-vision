%����POEģ����ţ��ŷ�������㷨�Ļ������涯��ѧ���
%Ŀǰֻ֧�����ֽ�
%Tor=idyna_poe_dh(robot,q,qd,qdd,gravity,f_end)
%����robotΪ������ģ�ͣ��������˶�ѧ�����Ͷ���ѧ����
%����qΪ�ؽڽǶȣ�1XN����qdΪ�ؽڽǼ��ٶ�,qddΪ��Ӧ�ؽڽǼ��ٶ�
%gravityΪ�������ٶȣ�f_endΪĩ��ִ����ʩ���غɣ�
%�������ڲ�����ʱ��ѡ��Ĭ��ֵ
%���TorqueΪ�ؽ�Ť��
%�ο�����Ϊ�������׵�.������ѧ��ģ�滮�����,7.8��:2018
%A Lie Group Formulation of Robot Dynamics

function Tor=idyna_poe_mdh(robot,q,qd,qdd,gravity,f_end)
%% ���úͻ�ȡ��������
%ȷ������Ƕȣ����ٶȣ��Ǽ��ٶ�Ϊ������
q=q(:)';
qd=qd(:)';
qdd=qdd(:)';
n=robot.n; %������������Ŀ
z0=[0,0,1]'; %������ϵ��z0�ķ���
grav= robot.gravity; %�������ٶȣ���������
fend = zeros(6, 1);%ĩ��ʩ�ӵ��غ�,W=[fx Fy Fz Mx My Mz]';
debug1=1;%���õ��Բ���
debug2=0;
if nargin>4%�ж��Ƿ������������ٶ�
    grav=gravity;
end
if nargin==6%�ж��ⲿ�Ƿ�ʩ���غ�
    fend=f_end;
end
%% �ж������Ƿ����Ҫ��
if numcols(q)~=n||numcols(qd)~=n||numcols(qd)~=n||numrows(q)~=1||...
        numrows(qd)~=1||numrows(qdd)~= 1||length(fend)~=6
    error('������������')
end
 %% ��� 
 offset=robot.offset;
 n=robot.n;
 V0=zeros(6,1);%���ó�ʼ�ٶ�
 Vd0=zeros(6,1);%���ó�ʼ���ٶȣ�����ο������о��������ʣ���Ϊ���ٶ�û�б��������У�
 
 %% ���T0�͹��Ծ���
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

%% �������
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

%% �������
for i=n:-1:1
    if i==n
        F{i}=Adjoint_v(inv(eye(4)))'*Fend+M{i}*Vd{i}-ad_v(V{i})'*M{i}*V{i};
    else
        F{i}=Adjoint_v(inv(T{i+1}))'*F{i+1}+M{i}*Vd{i}-ad_v(V{i})'*M{i}*V{i};
    end
    tor(i)=twist(i,:)*F{i};
end

%% ���Ծ�����
    function  M=m_matrix(m,pci,I)
        M(1:3,1:3)=m*eye(3);
        M(1:3,4:6)=-m*skew_poe(pci);
        M(4:6,1:3)=m*skew_poe(pci);
        M(4:6,4:6)=I-m*skew_poe(pci)*skew_poe(pci);
    end
end

