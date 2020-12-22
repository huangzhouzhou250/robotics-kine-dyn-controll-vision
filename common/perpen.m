%line=perpen(twist1,twist2) ����������������Ĺ�����
%����Ϊ������ʽ����w+p�����ߵķ����������һ�㣩
%������������ĸ�ʽѡ��

%������ 2020.9.1
function line=perpen(twist1,twist2)
%% ���봦�� 
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
    error('�����ʽ����')
end

%% �ж�����ֱ�ߵĹ�ϵ
%��������
v1=cross(r1,w1);
v2=cross(r2,w2);
if abs((w1'*v2+w2'*v1))>10^-5   %�����߲�����
    w3=cross(w1,w2);
    w3=w3/norm(w3);
    r12=r2-r1;
    r_temp=r12-r12'*w3*w3;
    r20=r1+r_temp;
    tw1=Twist('R',w1,r1);
    tw2=Twist('R',w2,r20);
    r3=twistcross(tw1,tw2);
    r3=r3(1:3);
elseif norm(cross(w1,w2))>10^-5 %�������ཻ
    w3=cross(w1,w2);
    w3=w3/norm(w3);
    tw1=Twist('R',w1,r1);
    tw2=Twist('R',w2,r2);
    r3=twistcross(tw1,tw2);
    r3=r3(1:3);
else   %������ƽ��
    error('������ƽ�У������߲�Ψһ')
end

%% ���
if out_flag==1
    line=Twist('R',w3,r3);
else
    line(1:3)=w3;
    line(4:6)=r3;
end
end 