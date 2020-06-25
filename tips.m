%% tips
%% 当动力学矩阵的系数太长无法输出时，可以储存为txt的格式，再复制。
fid=fopen('h11.txt','wt');
 fprintf(fid,'%s\n',char(H(1,1)))

 %% 计时
 t1=clock;
 for i=1:100
     dyn_lagr_p560(qn,qn,qn);
 end
t2=clock;

etime(t2,t1)  %  0.14s
%% 读写txt
filename = 'step.txt';
q=textread(filename);
%存储char
fid = fopen('d:/a.txt', 'wt');
fprintf(fid, '%s\n', char(f));
fclose(fid);

%写入数据时添加符号和格式
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

dlmwrite%储存方式多选
%% 
m=pi;
roundn(m,3)
vpa(m,2)

%% 输出
fprintf