%test
%函数测试文件


%% twistcross 测试 测试OK
%输入测试
% w1=[0 0 1]';r1=[0 0 0]';
% w2=[1 0 0]';r2=[1 0 1]';
% twist1=Twist('R',w1,r1);
% twist2=Twist('R',w2,r2);
% p=twistcross(twist1,twist2)
% twist1=[w1;r1];
% twist2=[w2;r2];
% p=twistcross(twist1,twist2)

%平行测试
% w1=[0 0 1]';r1=[0 0 0]';
% w2=[1 0 0]';r2=[1 0 1]';
% twist1=Twist('R',w1,r1);
% twist2=Twist('R',w1,r2);
% p=twistcross(twist1,twist2)

%不共面直线测试
% w1=[0 0 1]';r1=[0 0 0]';
% w2=[1 0 0]';r2=[1 1 1]';
% twist1=Twist('R',w1,r1);
% twist2=Twist('R',w2,r2);
% p=twistcross(twist1,twist2)


%% poe2dh(以UR5作为例子)
% 
% clc
% clear 
% 
% ur5 = xml2robot3d_21('UR5.xml');
% 
% twist=ur5.twist;
% T0=ur5.T0;
% for i=1:size(twist,1)
%     l{i}=Twist(twist(i,:));
%     w(:,i)=l{i}.w;
% end
% 
% % 获取相邻关节的公垂线
% x(:,1)=[1 0 0]';
% z(:,1)=[0 0 1]';
% r(:,1)=[0 0 0]';
% 
% for i=1:5
%     z(:,i+1)=w(:,i+1);
%     if norm(cross(w(:,i),w(:,i+1)))<10^-5
%         r1=l{i}.pole;
%         r2=l{i+1}.pole;
%         wtemp1=cross(w(:,i),r2-r1);
%         wtemp2=-cross(w(:,i),wtemp1);
%         x(:,i+1)=wtemp2/norm(wtemp2);
%         tw_temp=Twist('R',x(:,i+1),r(:,i));
%         r_temp=twistcross(tw_temp,l{i+1});
%         r(:,i+1)=r_temp(1:3);
%     else
%         per=perpen(l{i},l{i+1});
%         x(:,i+1)=per.w;
%         r_temp=twistcross(per,l{i+1});
%         r(:,i+1)=r_temp(1:3);
%     end
% end
% z(:,7)=T0(1:3,3);
% x(:,7)=T0(1:3,1);
% r(:,7)=T0(1:3,4);
% 
% for i=1:6
%     r12=r(:,i+1)-r(:,i);
%     d=r12'*z(:,i);
%     a=r12'*x(:,i+1);
%     alpha=acos(z(:,i)'*z(:,i+1));
%     offset=acos(x(:,i)'*x(:,i+1));
%     dh(i,1)=alpha;
%     dh(i,2)=a;
%     dh(i,3)=d;
%     dh(i,4)=offset;
%     L(i)=Link('alpha',alpha,'a',a,'d',d,'offset',offset);
% end
% 
% ur50=SerialLink(L)
% dh2=poe2dh(twist,T0)
%% 求两条线的公垂线 perpen(twist1,twist2)
% clc 
% clear
% 
% r1=[2,1,0];
% w1=[ 1 0 0];
% r2=[1 2 0];
% w2=[0 1 0];
% w3= cross(w1,w2);
% l1=Twist('R',w1,r1);
% l2=Twist('R',w2,r2);
% perpen(l1,l2)

% r12=r2-r1;
% r_temp=r12-r12*w3'*w3;
% r20=r1+r_temp;
% twist1=Twist('R',w1,r1);
% twist2=Twist('R',w2,r20);
% twistcross(twist1,twist2)


