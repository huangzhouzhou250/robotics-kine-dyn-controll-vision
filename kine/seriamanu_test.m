%seriamanu test
mdl_puma560;
p560p=dh2poe(p560);
q=zeros(1,6);

% 测量角度旋向是否正确
% for i=1:6
%     q(i)=1;
%     Tp=p560p.fkinep(q);
%     T=double(p560.fkine(q));
%     T-Tp
% end

%随机点,判断正解是否正确
% q1=rand(10,1)*4*pi-2*pi;
% q2=rand(10,1)*4*pi-2*pi;
% q3=rand(10,1)*4*pi-2*pi;
% q4=rand(10,1)*4*pi-2*pi;
% q5=rand(10,1)*4*pi-2*pi;
% q6=rand(10,1)*4*pi-2*pi;
% 
% for i=1:10
%     q=[q1(i) q2(i) q3(i) q4(i) q5(i) q6(i)];
%     Tp=p560p.fkinep(q);
%     T=double(p560.fkine(q));
%     T-Tp
% end

%