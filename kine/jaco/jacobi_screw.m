function J=jacobi_screw(twist,T0,q,argument)
%����������������ſ˱Ⱦ���
%��������twist����ʼλ��T0�͹ؽڽǶ�q
%����twistΪһ��6xN�ľ���NΪ�ؽ�����ǰ����������ת�ķ���
%������Ϊ������һ�㣻
%����T0Ϊĩ��ִ��������ʼλ�˾���qΪ�ؽڽǶȣ���λΪ����
%���Ϊ�ſ˱Ⱦ���J
%% �ж������Ƿ����Ҫ��
[twrow,twcol]=size(twist);
if (twrow~=6)||(twcol~=length(q))
    error('����Ĳ������淶')
end
%% ����������任����
J=zeros(6,twcol);
tw=cell(1,twcol);
T=cell(1,twcol);
for i=1:twcol
    tw{i}=Twist('R',twist(1:3,i),twist(4:6,i));
    T{i}=tw{i}.T(q(i));
    screw=tw{i}.double';
    if i>1 %���ʵʱ����
        T{i}=T{i-1}*T{i};
        screw=SE3(T{i-1}).Ad*screw; 
    end
    J(:,i)=screw;
end
if strcmp(argument,'s')
    return
elseif strcmp(argument,'b') %��������ſɱȾ���
    Tq=T{twcol}*T0;
    J=SE3(Tq).Ad*J;
end
end