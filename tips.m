%% tips
%% ������ѧ�����ϵ��̫���޷����ʱ�����Դ���Ϊtxt�ĸ�ʽ���ٸ��ơ�
fid=fopen('h11.txt','wt');
 fprintf(fid,'%s\n',char(H(1,1)))
 %% ��ʱ
 t1=clock;
 for i=1:100
     dyn_lagr_p560(qn,qn,qn);
 end
t2=clock;

etime(t2,t1)  %  0.14s
%% ��дtxt
filename = 'step.txt';
q=textread(filename);
%�洢char
fid = fopen('d:/a.txt', 'wt');
fprintf(fid, '%s\n', char(f));
fclose(fid);

%д������ʱ��ӷ��ź͸�ʽ
[m, n] = size(k);
fid=fopen('k.txt', 'wt');
for i=1:m
    for j=1:n
        if mod(j,3)~=0
            fprintf(fid,'%.3f,\t',k(i,j));
        else
            fprintf(fid,'%.3f\t',k(i,j));
        end
    end
   fprintf(fid,'\n');
end
fclose(fid);

dlmwrite%���淽ʽ��ѡ
%% 
m=pi;
roundn(m,3)
vpa(m,2)

%% ���
fprintf