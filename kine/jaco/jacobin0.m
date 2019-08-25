function J=jacobin0(DH,q)
%�������˵��ſɱȾ�������ΪDH����DH�͹ؽڽǶ�q
%DH����Ϊ��׼DH������ģ,��Ե���ȫ��������ת�ؽڵĻ�����
%
%% �ж������Ƿ����Ҫ��
[dhrow,dhcol]=size(DH);
if (dhcol~=4)||(dhrow~=length(q))
    error('����Ĳ���������Ҫ��')
end
%���ӻ����굽�ؽ�i����ı任����
T=cell(1,6);%����Ԫ������
for i=1:dhrow
    T{i}=DHstd(DH(i,1),DH(i,2),DH(i,3),DH(i,4)+q(i),'rad');%����������˱任����
    if i>1
       T{i}=T{i-1}*T{i}; %��ǰһ�������ȡ0��i�ı任����
    end
end
%% ���ſ˱Ⱦ���
%Ԥ����ȡz0,pe
z0=[0 0 1]';
pe=T{6}*[0 0 0 1]';
pe=pe(1:3);
J=zeros(6,dhrow);

for i=1:dhrow
    %���pi��zi
    if i==1
        w=z0;
        p=[ 0 0 0]';
    else
        w=t2r(T{i-1})*z0;
        p=T{i-1}*[0 0 0 1]';
        p=p(1:3);
    end
    v=cross(w,pe-p);
    %�ϳ��ſ˱Ⱦ���
    J(:,i)=[v;w];
end
end