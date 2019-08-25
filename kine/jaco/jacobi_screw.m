function J=jacobi_screw(twist,T0,q,argument)
%利用旋量方法求解雅克比矩阵
%利用旋量twist，起始位置T0和关节角度q
%输入twist为一个6xN的矩阵，N为关节数，前三行轴线旋转的方向，
%后三行为轴线上一点；
%输入T0为末端执行器的起始位姿矩阵，q为关节角度，单位为弧度
%输出为雅克比矩阵J
%% 判断输入是否符合要求
[twrow,twcol]=size(twist);
if (twrow~=6)||(twcol~=length(q))
    error('输入的参数不规范')
end
%% 利用旋量求变换矩阵
J=zeros(6,twcol);
tw=cell(1,twcol);
T=cell(1,twcol);
for i=1:twcol
    tw{i}=Twist('R',twist(1:3,i),twist(4:6,i));
    T{i}=tw{i}.T(q(i));
    screw=tw{i}.double';
    if i>1 %求解实时旋量
        T{i}=T{i-1}*T{i};
        screw=SE3(T{i-1}).Ad*screw; 
    end
    J(:,i)=screw;
end
if strcmp(argument,'s')
    return
elseif strcmp(argument,'b') %求解物体雅可比矩阵
    Tq=T{twcol}*T0;
    J=SE3(Tq).Ad*J;
end
end