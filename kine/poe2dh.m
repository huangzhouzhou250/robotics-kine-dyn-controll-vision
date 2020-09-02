%poe2dh��twist,T0��
%����֪������е�۵�ָ����ģ��ת��Ϊ��׼��DHģ��
%twistΪ�����������Ԫ�����飬���߻��ھ������ʽ���
%ÿ�б�ʾһ������������ǰ����Ϊv������Ϊw

%�ο����ף���������ѧ����ģ���滮����ơ�
%��������ѧ����ģ�滮����ơ�

%������ 2020.9.2

function dh=poe2dh(twist,T0)
%% �жϲ���������
if isa(twist,'double') && size(twist,2)==6   %���������ʽ����
    n=size(twist,1);
    for i=1:n
        l{i}=Twist(twist(i,:));
    end
elseif isa(l,'cell') && isa(l{1},'Twist')    %������Ԫ���������ʽ����
    n=length(twist);
    l=twist;
end

%% ȷ��ÿ���ؽ�����ϵ
%����ϵ1
x(:,1)=[1 0 0]'; 
z(:,1)=[0 0 1]';
r(:,1)=[0 0 0]';    %����ϵԭ��

for i=1:5
    z_temp=l{i+1}.w;        %��ȡ������w�����ж��Ƿ�Ϊ��ת����
    if norm(z_temp)>10^-5   %��ת��
        z(:,i+1)=z_temp;    %����ת�᷽����Ϊ����ϵz��
        if norm(cross(l{i}.w,l{i+1}.w))<10^-5   %�ж��Ƿ�����ƽ�У���������ƽ��ʱ�������߲�Ψһ
            r1=l{i}.pole;   %���������һ��
            r2=l{i+1}.pole;
            wtemp1=cross(l{i}.w,r2-r1);
            wtemp2=-cross(l{i}.w,wtemp1);   
            x(:,i+1)=wtemp2/norm(wtemp2);   %ͨ���������ȷ��ƽ�����ߵĹ����߷���
            tw_temp=Twist('R',x(:,i+1),r(:,i));  %ʹ�����߾�����һ����ϵ��ԭ�㣬�򻯼���
            r_temp=twistcross(tw_temp,l{i+1});   %�������ϵԭ��
            r(:,i+1)=r_temp(1:3);
        else                %���߲�ƽ��ʱ
            per=perpen(l{i},l{i+1});    %��⹫���ߣ�����μ�����perpen
            x(:,i+1)=per.w;
            r_temp=twistcross(per,l{i+1});  %��⽻��
            r(:,i+1)=r_temp(1:3);
        end
    else                    %�ƶ��ؽ�
        z(:,i+1)=l{i+1}.v;  %���߷�����Ϊz��
        r(:,i+1)=r(:,i);    %����һ����ϵ��ԭ�㶨��Ϊ�ƶ����ߵ�ԭ��
        x_temp=cross(z(:,i),z(:,i+1));   
        x=x_temp/norm(x_temp);   %ȷ��x�᷽��
    end
end
%���ĩ������ϵ
z(:,n+1)=T0(1:3,3);
x(:,n+1)=T0(1:3,1);
r(:,n+1)=T0(1:3,4);

%% ���DH������ 
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
    if  norm(l{i}.w)>10^-5  %���ڱ���Ƿ�Ϊ��ת�ؽ�
        dh(i,5)=1;
    else
        dhdh(i,5)=0;
    end
end

end