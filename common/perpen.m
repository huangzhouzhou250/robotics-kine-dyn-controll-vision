%line=perpen(twist1,twist2) 用于求解两个旋量的公垂线
%输入为旋量格式或者w+p（轴线的方向和轴线上一点）
%输出会根据输入的格式选择

%黄洲洲 2020.9.1
function line=perpen(twist1,twist2)
%% 输入处理 
if isa(twist1,'Twist') && isa(twist2,'Twist')
    w1=twist1.w;
    r1=twist1.pole;
    w2=twist2.w;
    r2=twist2.pole;
    w1=w1(:);
    r1=r1(:);
    w2=w2(:);
    r2=r2(:);
    out_flag=1;
elseif isa(twist1,'double') && isa(twist2,'double') && length(twist1)==6 ...
        && length(twist2)==6
    w1=twist1(1:3);
    r1=twist1(4:6);
    w2=twist2(1:3);
    r2=twist2(4:6);
    w1=w1(:);
    r1=r1(:);
    w2=w2(:);
    r2=r2(:);
else
    error('输入格式不对')
end

%% 判断两条直线的关系
%旋量互矩
v1=cross(r1,w1);
v2=cross(r2,w2);
if abs((w1'*v2+w2'*v1))>10^-5   %两轴线不共面
    w3=cross(w1,w2);
    w3=w3/norm(w3);
    r12=r2-r1;
    r_temp=r12-r12'*w3*w3;
    r20=r1+r_temp;
    tw1=Twist('R',w1,r1);
    tw2=Twist('R',w2,r20);
    r3=twistcross(tw1,tw2);
    r3=r3(1:3);
elseif norm(cross(w1,w2))>10^-5 %两轴线相交
    w3=cross(w1,w2);
    w3=w3/norm(w3);
    tw1=Twist('R',w1,r1);
    tw2=Twist('R',w2,r2);
    r3=twistcross(tw1,tw2);
    r3=r3(1:3);
else   %两轴线平行
    error('两条线平行，公垂线不唯一')
end

%% 输出
if out_flag==1
    line=Twist('R',w3,r3);
else
    line(1:3)=w3;
    line(4:6)=r3;
end
end 