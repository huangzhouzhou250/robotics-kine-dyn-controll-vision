%poe2dh（twist,T0）
%将已知串联机械臂的指数积模型转换为标准的DH模型
%twist为旋量数组或者元胞数组，或者基于矩阵的形式表达
%每行表示一个旋量，其中前三列为v后三行为w

%参考文献：《机器人学：建模、规划与控制》
%《机器人学：建模规划与控制》

%黄洲洲 2020.9.2

function dh=poe2dh(twist,T0)
%% 判断并调整输入
if isa(twist,'double') && size(twist,2)==6   %以数组的形式输入
    n=size(twist,1);
    for i=1:n
        l{i}=Twist(twist(i,:));
    end
elseif isa(l,'cell') && isa(l{1},'Twist')    %以旋量元胞数组的形式输入
    n=length(twist);
    l=twist;
end

%% 确定每个关节坐标系
%坐标系1
x(:,1)=[1 0 0]'; 
z(:,1)=[0 0 1]';
r(:,1)=[0 0 0]';    %坐标系原点

for i=1:5
    z_temp=l{i+1}.w;        %获取旋量的w用于判断是否为旋转轴线
    if norm(z_temp)>10^-5   %旋转轴
        z(:,i+1)=z_temp;    %将旋转轴方向定义为坐标系z轴
        if norm(cross(l{i}.w,l{i+1}.w))<10^-5   %判断是否相邻平行，相邻轴线平行时，公垂线不唯一
            r1=l{i}.pole;   %求解轴线上一点
            r2=l{i+1}.pole;
            wtemp1=cross(l{i}.w,r2-r1);
            wtemp2=-cross(l{i}.w,wtemp1);   
            x(:,i+1)=wtemp2/norm(wtemp2);   %通过两个叉积确定平行轴线的公垂线方向
            tw_temp=Twist('R',x(:,i+1),r(:,i));  %使公垂线经过上一坐标系的原点，简化计算
            r_temp=twistcross(tw_temp,l{i+1});   %求解坐标系原点
            r(:,i+1)=r_temp(1:3);
        else                %轴线不平行时
            per=perpen(l{i},l{i+1});    %求解公垂线，具体参见函数perpen
            x(:,i+1)=per.w;
            r_temp=twistcross(per,l{i+1});  %求解交点
            r(:,i+1)=r_temp(1:3);
        end
    else                    %移动关节
        z(:,i+1)=l{i+1}.v;  %轴线方向定义为z轴
        r(:,i+1)=r(:,i);    %将上一坐标系的原点定义为移动轴线的原点
        x_temp=cross(z(:,i),z(:,i+1));   
        x=x_temp/norm(x_temp);   %确定x轴方向
    end
end
%求解末端坐标系
z(:,n+1)=T0(1:3,3);
x(:,n+1)=T0(1:3,1);
r(:,n+1)=T0(1:3,4);

%% 求解DH参数表 
for i=1:6
    r12=r(:,i+1)-r(:,i); 
    d=r12'*z(:,i);   %
    a=r12'*x(:,i+1);
    alpha=acos(z(:,i)'*z(:,i+1));
    theta=acos(x(:,i)'*x(:,i+1));
    dh(i,1)=alpha;
    dh(i,2)=a;
    dh(i,3)=d;
    dh(i,4)=theta;
    if  norm(l{i}.w)>10^-5  %用于标记是否为旋转关节
        dh(i,5)=1;
    else
        dhdh(i,5)=0;
    end
end

end