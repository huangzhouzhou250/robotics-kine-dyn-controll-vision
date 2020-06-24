%�˺���������ⷴ�Գƾ���ͽ�����ת��Ϊ4X4����
%twist ����Ϊ1X3����1X6������
%���seΪ3X3����4X4�ľ���
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
    error('����������������')
end
end