%��������ָ������ʽ�Ĵ���������ģ��
%��������Ӧ�ò��������˺����λ����˵Ĵ���
%����ָ���������ݿ��Բο�һ�����ף�
%�����׵ȡ�������ѧ����ģ�������Ӿ���
%������ȡ�������ѧ�ļ��λ�����
%������Բο�Peter Corke��RTB������

classdef SerialManu < handle 
     properties
         name       %����������
         comment    %���������˵��
         T0         %�趨������ĩ������ڻ�����ϵ�ĳ�ʼλ��
         gravity    %�������ٶ�ʸ��
         base       %�����˻�������ϵ�任
         tool       %�����˹�������ϵ�任
         plot3dopt     %3d����
         ikineType  %�������
         faces      %��ȡ�����ļ�֮�����͵�
         points
     end
     events
         moved
     end
     properties (SetAccess = private)
         n
         jolinks
         T
     end
     properties (Dependent = true)
         %�˴�Ϊ�������ԣ��������Բ���ͨ����ֵ����
         %ֻ����get��������
         w
         r
         v
         twist
         offset
         qlim
         qdlim
         qddlim
         torlim
         theta
         serialtype  %����������,��'RRRRRR'

     end
     
     methods
         function robot=SerialManu(varargin)
             %����������е�۶���
             % r=SerialManu()���ڴ���һ���յĶ���
             %���и������ѡȡΪĬ��ֵ
             %
             %r=SerialManu(robot)���ڸ���һ�������˶���
             %����ȫ������������
             %
             %r=SerialManu([jl1 jl2 ... ]��option)���ò�������������е��
             %[jl1 jl2 ... ]��ʾ�ؽ��������ΰڷ�
             %
             %options:
             %'name',name            ����������
             %'comment',comment      �����˵����˵������������˳��̣�ʱ�䣬���Եȵ�
             %'base',T               �趨�����˻����ı任����4X4double or SE3��
             %'tool',T               �趨�����˵Ĺ�������ϵ�任��4X4double or SE3��
             %'gravity',G            ��������,����Ϊ��ά����
             %plot3dopt,p               ������άģ�͵ķ�ʽ��n�޷����ƣ�y���Ի��ƣ�
             
             %����Ĭ�ϲ���
             opt.jolinks=[];
             opt.n=0;
             opt.name = '';
             opt.comment ='';
             opt.base = eye(4);
             opt.tool = eye(4);
             opt.gravity = [0; 0; 9.81];
             opt.T0=[];
             opt.plot3dopt = {};
             opt.ikine = {};

             
             [opt,arg] = tb_optparse(opt, varargin);
                 
             if isempty(varargin)
                 %��Ĭ�ϲ�������Ĭ�϶���;
                 robot.n=0;
                 robot.jolinks=[];
                 robot.name='';
                 robot.comment = '';
                 robot.T0=eye(4);
                 robot.base = eye(4);
                 robot.tool=eye(4);
                 robot.gravity=[0; 0; 9.81];
                 robot.plot3dopt={};
             elseif nargin==1 && isa(varargin{1},'SerialManu')
                 %���ƴ��������˶���
                 this = varargin{1};
                 robot.n=this.n;
                 robot.jolinks=this.jolinks;
                 robot.name=this.name;
                 robot.comment = this.comment;
                 robot.T0=this.T0;
                 robot.base = this.base;
                 robot.tool=this.tool;
                 robot.gravity=this.gravity;       
                 robot.plot3dopt=this.gravity;
             elseif nargin>=2 &&  isa(arg{1},'JoLink')
                 %���ݲ�����������������
                 Link=arg{1};
                 robot.jolinks=Link;
                 robot.n=length(robot.jolinks);
                 robot.name=opt.name;
                 robot.comment = opt.comment;
                 robot.T0=opt.T0;
                 robot.base = opt.base;
                 robot.tool=opt.tool;
                 robot.gravity=opt.gravity;
             else
                 error('�����ʽ����');
             end
             
         end%���캯��
         
         
         function display(robot)
             %��������˲���
             
         end
         %Ϊ�������Զ���set����
         function k=get.w(robot)
             if robot.n==0
                 k=[];
             else
                 
                 k = [robot.jolinks.w]';
             end
         end %w����
         
         function k=get.v(robot)
             if robot.n==0
                 k=[];
             else
                 k = [robot.jolinks.v]';
             end
         end %v����
         
         function k=get.offset(robot)
             if robot.n==0
                 k=[];
             else
                 k = [robot.jolinks.offset];
             end
         end %offset����
         
         function k=get.qlim(robot)
             if robot.n==0
                 k=[];
             else
                 for i=1:robot.n
                     k(i,:)=robot.jolinks(i).qlim;
                 end
             end
         end %qlim����
         
         function k=get.qdlim(robot)
             if robot.n==0
                 k=[];
             else
                 k = [robot.jolinks.qdlim];
             end
         end %qdlim����
         
         function k=get.qddlim(robot)
             if robot.n==0
                 k=[];
             else
                 k = [robot.jolinks.qddlim];
             end
         end %qddlim����
         
         function k=get.torlim(robot)
             if robot.n==0
                 k=[];
             else
                 k = [robot.jolinks.torlim];
             end
         end %w����
         
         function k=get.twist(robot)
             k=[robot.w robot.v];
         end %w����
         
         function k=get.serialtype(robot)
             if robot.n==0
                 k=[];
             else
                 k = [robot.jolinks.jointtype];
             end
         end %w����
         
         
         
     end%methods
end%��