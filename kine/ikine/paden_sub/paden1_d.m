function [ theta] = paden1_d( twist,p,q,De)
%Paden3 ����������һ�㾭����ת����ת�����q����ΪDe��λ��
%twistΪ��ת�ؽ��������������pΪ���
%theta1��thetaΪ���ܵĽǶ�
%% ������
r=twist.pole;
w=twist.w;
u=p-r;
v=q-r;
u1=u-w*w'*u;
v1=v-w*w'*v;
De_0=De^2-(w'*(p-q))^2;
De1=sqrt(De_0);
theta_0=atan2(w'*cross(u1,v1),u1'*v1);
t=acos((norm(u1)^2+norm(v1)^2-De1^2)/(2*norm(u1)*norm(v1)));
q1=theta_0+t;
q=theta_0-t;
theta=[q q1];
% ����Ƕȳ�������
for i=1:2
    if q(i)>pi
        q(i)=q(i)-2*pi;
    elseif q(i)<-pi
        q(i)=q(i)+2*pi;
    end
end

end
