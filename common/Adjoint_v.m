%
function T6=Adjoint_v(T4)
if size(T4,1)~=4 || size(T4,2)~=4 
    error('输入参数有误')
end
T6=zeros(6,6);

R=T4(1:3,1:3);
p=T4(1:3,4);

T6(1:3,1:3)=R;
T6(4:6,4:6)=R;
T6(1:3,4:6)=skew_poe(p)*R;
end