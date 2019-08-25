function J=jacobin0(DH,q)
%求解机器人的雅可比矩阵，输入为DH参数DH和关节角度q
%DH参数为标准DH方法建模,针对的是全部都是旋转关节的机器人
%
%% 判断输入是否符合要求
[dhrow,dhcol]=size(DH);
if (dhcol~=4)||(dhrow~=length(q))
    error('输入的参数不符合要求')
end
%求解从基坐标到关节i坐标的变换矩阵
T=cell(1,6);%声明元胞数组
for i=1:dhrow
    T{i}=DHstd(DH(i,1),DH(i,2),DH(i,3),DH(i,4)+q(i),'rad');%求解相邻连杆变换矩阵
    if i>1
       T{i}=T{i-1}*T{i}; %与前一个相乘求取0到i的变换矩阵
    end
end
%% 求雅克比矩阵
%预先求取z0,pe
z0=[0 0 1]';
pe=T{6}*[0 0 0 1]';
pe=pe(1:3);
J=zeros(6,dhrow);

for i=1:dhrow
    %求解pi和zi
    if i==1
        w=z0;
        p=[ 0 0 0]';
    else
        w=t2r(T{i-1})*z0;
        p=T{i-1}*[0 0 0 1]';
        p=p(1:3);
    end
    v=cross(w,pe-p);
    %合成雅克比矩阵
    J(:,i)=[v;w];
end
end