function Tq=fkine_screw(twist,T0,q)
%��������twist����ʼλ��T0�͹ؽڽǶ�q
%����twistΪһ��6xN�ľ���NΪ�ؽ�����ǰ����������ת�ķ���
%������Ϊ������һ�㣻
%����T0Ϊĩ��ִ��������ʼλ�˾���qΪ�ؽڽǶȣ���λΪ����
%���TqΪ��������תq��ĩ��ִ������λ��
%% �ж������Ƿ����Ҫ��
[twrow,twcol]=size(twist);
if (twrow~=6)||(twcol~=length(q))
    error('����Ĳ������淶')
end
%% ����������任����
tw=cell(1,twcol);
T=cell(1,twcol);
for i=1:twcol
    tw{i}=Twist('R',twist(1:3,i),twist(4:6,i));
    T{i}=tw{i}.T(q(i));
    if i>1
        T{i}=T{i-1}*T{i};
    end
end
%% ��Tq
Tq=T{twcol}*T0;
end