%指数积模型中连杆的建立
%关于指数积公式的相关内容参考文献如下：
%熊有伦等《机器人学：建模控制与视觉》
%李泽湘等《机器人学的几何基础》
%程序参考了Peter Corke的RTB工具箱


classdef JoLink < matlab.mixin.Copyable
    properties
        %运动学参数
        w   %关节轴线位置
        r   %关节轴线上一点
        jointtype  %关节类型，revolute='R', prismatic='P' -- should be an enum
        offset%关节默认偏角
        name  %连杆名字
        flip %关节反向移动
        qlim %关节角度限制(2x1)
        qdlim %关节角速度限制
        qddlim %关节角加速度限制
        % 动力学参数
        m  % l连杆质量
        rc % 连杆质心相对于基坐标系的位置
        I  % 连杆相对于自身质心的惯性矩阵
        %电机参数
        Jm % 电机等效转动惯量
        B  % 电机粘性摩擦系数 (1x1 or 2x1)
        torlim %关节力矩限制
        Tc % 电机库伦摩擦系数 (1x2 or 2x1)
        G  % 电机减速比  
    end %类的性质
    methods
        function jl=JoLink(varargin)
            %JoLink的构造函数，用于创建连杆对象
            %
            %创建方式如下：
            %jl=JoLink();按照默认参数创建连杆
            %jl=JoLink(lolink);复制其他连杆对象
            %jl=JoLink(options);根据运动学和动力学参数创建连杆
            %options::
            %'w',w;关节轴线在基坐标系下的方向向量
            %'r',r;关节轴线上一点在基坐标系中的位置
            %'offset',o;默认起始关节角，用于调整使得模型与实物角度一致(defult 0)
            %'flip',logic;为true时关节旋向相反，false一致(default false)
            %'qlim',[a,b];关节角度限制（1X2,default[-pi,pi]）
            %'qdlim',qd;关节角速度限制（default[]）
            %'qddlim',qd;关节角加速度限制（default[]）
            %'m'，m;连杆质量（default 0）
            %'rc',rc;连杆质心相对于基坐标系的位置向量（default []）
            %'I',i;连杆相对于其质心的惯性矩阵(default [])
            % 'Jm',J  ;电机转动惯量(default 0)
            % 'B',B  ;关节粘滞摩擦系数（default 0）
            %'torlim',tor %关节力矩限制
            % 'Tc',T  ;电机库仑摩擦系数 (1x1 or 2x1), (default [0 0])
            % 'G',G ；电机减速比(default 101)
            if nargin==0
                %创建默认连杆
                
                %运动学参数
                jl.w=[0 0 1]';
                jl.r=[0 0 0]';
                jl.offset=0;
                jl.flip=false;
                jl.name='';
                jl.qlim=[-pi pi];
                jl.qdlim=pi/2;%rad/s
                jl.qddlim=0.5;%rad/s^2
                jl.jointtype='R';
                %动力学参数
                jl.m=0;
                jl.rc=[];
                jl.I=[];
                %电机参数
                jl.Jm=0;
                jl.B=0;
                jl.torlim=0;
                jl.Tc=0;
                jl.G=101;
            elseif nargin==1 && isa(varargin{1},'JoLink')  %需要提取出元胞数组中的连杆
                %复制对象
                this = varargin{1};
                for j=1:length(this)
                    jl(j) = JoLink();
                    %复制属性
                    p = properties(this(j));
                    for i = 1:length(p)
                        jl(j).(p{i}) = this(j).(p{i});
                    end
                end
            else
                %根据给定参数创建连杆
                %设定相关参数
                opt.w= [];
                opt.r = [];
                opt.G = 0;
                opt.B = 0;
                opt.Tc = [0 0];
                opt.Jm = 0;
                opt.I = [];
                opt.m = 0;
                opt.rc = [0 0 0];
                opt.offset = 0;
                opt.qlim = [-pi pi];
                opt.qdlim =-pi/2;
                opt.qddlim =0.5;
                opt.type = {'revolute', 'prismatic', 'fixed'};
                opt.flip = false;
                opt.torlim=0;
                
                [opt,args] = tb_optparse(opt, varargin);
                
                if isempty(args)
                    %判断连杆类型
                    if norm(opt.w)==0
                        jl.jointtype='P';
                    else
                        jl.jointtype='R';
                    end
                    jl.w= opt.w;
                    jl.r= opt.r;
                    jl.offset =opt.offset;
                    jl.flip =opt.flip;
                    
                    jl.qlim =  opt.qlim;
                    jl.qdlim =  opt.qdlim;
                    jl.qddlim =   opt.qddlim;
                    
                    jl.m =opt.m;
                    jl.r = opt.r;
                    jl.I = opt.I;
                    jl.Jm =opt.Jm;
                    jl.G =opt.G;
                    jl.B =opt.B;
                    jl.Tc = opt.Tc;
                    jl.torlim=opt.torlim;
                end
            end
        end%构造函数
        
        function tor=friction(jl,qd)
            %关节摩擦力求解
            %参考Link中同名函数编写
            tor = jl.B * abs(jl.G) * qd;
            % Coulomb friction
            if ~isa(qd, 'sym')
                if qd > 0
                    tor = tor + jl.Tc(1);
                elseif qd < 0
                    tor = tor + jl.Tc(2);
                end
            end
            % scale up by gear ratio
            tor = -abs(jl.G) * tor;     % friction opposes motion
        end % 摩擦力
        
        function logic=islimt(jl,q,qd,qdd,tor)
            %用于判断关节角度，关节角速度，角加速度和力是否超出限制
        end
    end %method
end %类的结束