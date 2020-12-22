function [ theta] = paden1_d( twist,p,q,De)
%Paden3 问题描述：一点经过旋转轴旋转到离点q距离为De的位置
%twist为旋转关节所代表的旋量，p为起点
%theta1和theta为可能的角度
%% 求解过程
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
% 避免角度超出限制
for i=1:2
    if q(i)>pi
        q(i)=q(i)-2*pi;
    elseif q(i)<-pi
        q(i)=q(i)+2*pi;
    end
end

end
