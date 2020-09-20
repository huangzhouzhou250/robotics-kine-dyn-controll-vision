% ikine_num
function [q,k,dis]=ikine_num(robot,Tg,q0)
%�������˵���ֵ��,ֻ�����ݵ����������һ��⣻
%����޷�������н⣬��ʱ��Ҳ�޷�����Լ���Ҫ�Ľ⡣
%��Ϊ���԰汾���ó����а���������Ĳ�����
%����ΪrobotΪ�����˽ṹ�壬���������˵��˶�ѧ������
%TgΪ����ξ������Ŀ��λ�ˣ�Ϊ4X4double����SE3
%q0Ϊ��������ʼ�㣬��������룬��Ϊ[0 0 0 0 0 0];
n=robot.n; %��ȡ������������
%�趨��ʼ�Ƕ�
qi=zeros(1,n);
if nargin==3
    qi=q0;
end
if length(qi)~=n
    error('����Ƕ�����')
end
%ת��TgΪ4X4double��������ֲ����������
if isa(Tg,'SE3')
    Tg=Tg.T;
end
k=1;
flag=true;%����ѭ����־
dis2=1;
while flag
    T0=robot.fkine(qi);%����ʼ״̬ʱλ�˾���
    T0=T0.T;
    T=Tg/T0;
    v=T(1:3,4);         %��������任��Ӧ��ƽ����
    w=vex(t2r(T)-eye(3)); %����������ת����ת��Ϊ��ǵĳ˻�
    delta=[v;w];        %�ϲ�Ϊ6άʸ��
    dis1=norm(delta);    %��deltaģ����Ϊ�ж�ֵ
    jaco0=robot.jacob0(qi);  %��������q0ʱ���ſ˱Ⱦ���
    if n<6
        jaco=jaco0'*jaco0;
    else
        jaco=jaco0*jaco0';
    end
    det_j=sqrt(abs(det(jaco)));   %����ſ˱�����ʽ�ľ���ֵ
     %����det_j��ȡ��ͬ�ĵ�������
    if  det_j>10^-4
        if n<6
            jaco_inv=inv(jaco)*jaco0';
        else  
            jaco_inv=jaco0'*inv(jaco);
        end
        temp=jaco_inv*delta;
        alpha=1;
    else
        %����λ�ø������õķ���
        lambda=0.01;
        temp=jaco'*inv(jaco+lambda*eye(n))*delta;
        alpha=1.5;
    end
    if dis1>0.1
        beta=0.3;
    else
        beta=0.92;
    end
    
    qi=qi+alpha*temp'*beta;          %�����º�ĽǶ�
    while ~(min(qi-ones(1,n)*-pi)>0 && min(ones(1,n)*pi-qi)>0)
        for i=1:length(qi)
            if qi(i)>=pi
                qi(i)=qi(i)-2*pi;
            elseif qi(i)<-pi
                qi(i)=qi(i)+2*pi;
            end
        end
    end
    if ((abs(dis2-dis1))<10^-6)||(k>200)        %�ж��Ƿ��������
        flag=false;
    end
    dis2=dis1;
    if nargout>1
    dis(k)=dis2;
    k=k+1;
    end
end
q=qi;
% plot(dis)
end