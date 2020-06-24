function T=adjoint_v(twist)
twist=twist(:);
if length(twist)~=6
    error('输入旋量长度不对')
end
T=zeros(6,6);

w=twist(4:6);
v=twist(1:3);

w_skew=skew_poe(w);
v_skew=skew_poe(v);

T(1:3,1:3)=w_skew;
T(4:6,4:6)=w_skew;
T(1:3,4:6)=v_skew;
end