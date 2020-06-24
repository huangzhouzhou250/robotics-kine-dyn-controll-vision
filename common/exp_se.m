% 此函数是罗德里格斯公式，用于求解se3与SE3之间的指数关系
%twist为输入的旋量，[v,w]为六维矢量
%q为关节变量
%输出T为变换矩阵
function T=exp_se(twist,q)
T=eye(4);
twist=twist(:);
if length(twist)~=6
    error('输入旋量长度有误')
end
%% 获取w 和 v 
w=twist(4:6);
v=twist(1:3);
%% 求解
if abs(norm(w)-0)>10^-4   %旋转关节
    %归一化
    w=w/norm(w);
    v=v/norm(w);
    w_skew=skew_poe(w);
    R=eye(3)+w_skew*sin(q)+w_skew*w_skew*(1-cos(q));
    p=(eye(3)-R)*(w_skew*v)+w*w'*v*q;
    T(1:3,1:3)=R;
    T(1:3,4)=p;
elseif abs(norm(w)-0)<10^-4 && abs(norm(v)-0)>10^-4 %移动关节
    v=v/norm(v);
    p=v*q;
    T(1:3,4)=p;
else
    error('输入旋量格式不对')
end
end