% �˺������޵����˹��ʽ���������se3��SE3֮���ָ����ϵ
%twistΪ�����������[v,w]Ϊ��άʸ��
%qΪ�ؽڱ���
%���TΪ�任����
function T=exp_se(twist,q)
T=eye(4);
twist=twist(:);
if length(twist)~=6
    error('����������������')
end
%% ��ȡw �� v 
w=twist(4:6);
v=twist(1:3);
%% ���
if abs(norm(w)-0)>10^-4   %��ת�ؽ�
    %��һ��
    w=w/norm(w);
    v=v/norm(w);
    w_skew=skew_poe(w);
    R=eye(3)+w_skew*sin(q)+w_skew*w_skew*(1-cos(q));
    p=(eye(3)-R)*(w_skew*v)+w*w'*v*q;
    T(1:3,1:3)=R;
    T(1:3,4)=p;
elseif abs(norm(w)-0)<10^-4 && abs(norm(v)-0)>10^-4 %�ƶ��ؽ�
    v=v/norm(v);
    p=v*q;
    T(1:3,4)=p;
else
    error('����������ʽ����')
end
end