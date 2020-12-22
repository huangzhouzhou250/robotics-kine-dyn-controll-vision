function [theta1,theta2] = paden2in_d(twist1,twist2,p,q1,q2,de1,de2)
%Padene ������������p����������ཻ��ת�ᣨtwist2,twist1Ϊ������ߵ��˶�������
%����q1����de1��q2����Ϊde2�ĵ�
%r0Ϊ���������Ľ���
%���Ϊ������ܵĹؽڽǶ�
%% ������
%���ݴ���
r0=twistcross(twist1,twist2);
r0=r0(1:3);
p=p(:);
q1=q1(:);
q2=q2(:);
%���
w=p-r0;
u=q1-r0;
v=q2-r0;

de3=norm(w);
de4=norm(u);
de5=norm(v);

u1=u/(norm(u));
v1=v/(norm(v));

t1=de3*((de3^2+de4^2-de1^2)/(2*de3*de4));
t2=de3*((de3^2+de5^2-de2^2)/(2*de3*de5));

x1=(t1-u1'*v1*t2)/(1-(u1'*v1)^2);
x2=(t2-u1'*v1*t1)/(1-(u1'*v1)^2);
x3_1=(de3^2-x1^2-x2^2-2*x1*x2*u1'*v1)/(norm(cross(u1,v1))^2);
if abs(x3_1-0)<0.01*min(abs(x1),abs(x2))
    x3=0;
elseif x3_1>0
   x3=sqrt(x3_1);
else
    error('there is no answer');
end     
z=x1*u1+x2*v1+x3*cross(u1,v1);
q=z+r0;
[th1,th2]=paden2in(twist1,twist2,p,q);
z1=x1*u1+x2*v1-x3*cross(u1,v1);
q1=z1+r0;
[th3,th4]=paden2in(twist1,twist2,p,q1);
theta1=[th1 th3 ];
theta2=[th2 th4 ];
end
