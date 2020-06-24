%此函数用于求解反对称矩阵和将旋量转换为4X4矩阵
%twist 输入为1X3或者1X6的向量
%输出se为3X3或者4X4的矩阵
function se=skew_poe(twist)
twist=twist(:);
n=length(twist);
if n==3
    se=zeros(3,3);
    se(3,2)=twist(1);
    se(2,3)=-twist(1);
    se(3,1)=-twist(2);
    se(1,3)=twist(2);
    se(1,2)=-twist(3);
    se(2,1)=twist(3);
elseif n==6
    se=zeros(4,4);
    se(3,2)=twist(1);
    se(2,3)=-twist(1);
    se(3,1)=-twist(2);
    se(1,3)=twist(2);
    se(1,2)=-twist(3);
    se(2,1)=twist(3);
    se(1:3,4)=twist(4:6);
else
    error('输入向量长度有误')
end
end