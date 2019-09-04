%ָ����ģ�������˵Ľ���
%����ָ������ʽ��������ݲο��������£�
%�����׵ȡ�������ѧ����ģ�������Ӿ���
%������ȡ�������ѧ�ļ��λ�����
%����ο���Peter Corke��RTB������


classdef JoLink < matlab.mixin.Copyable
    properties
        %�˶�ѧ����
        w   %�ؽ�����λ��
        r   %�ؽ�������һ��
        jointtype  %�ؽ����ͣ�revolute='R', prismatic='P' -- should be an enum
        offset%�ؽ�Ĭ��ƫ��
        name  %��������
        flip %�ؽڷ����ƶ�
        qlim %�ؽڽǶ�����(2x1)
        qdlim %�ؽڽ��ٶ�����
        qddlim %�ؽڽǼ��ٶ�����
        % ����ѧ����
        m  % l��������
        rc % ������������ڻ�����ϵ��λ��
        I  % ����������������ĵĹ��Ծ���
        %�������
        Jm % �����Чת������
        B  % ���ճ��Ħ��ϵ�� (1x1 or 2x1)
        torlim %�ؽ���������
        Tc % �������Ħ��ϵ�� (1x2 or 2x1)
        G  % ������ٱ�  
    end %�������
    methods
        function jl=JoLink(varargin)
            %JoLink�Ĺ��캯�������ڴ������˶���
            %
            %������ʽ���£�
            %jl=JoLink();����Ĭ�ϲ�����������
            %jl=JoLink(lolink);�����������˶���
            %jl=JoLink(options);�����˶�ѧ�Ͷ���ѧ������������
            %options::
            %'w',w;�ؽ������ڻ�����ϵ�µķ�������
            %'r',r;�ؽ�������һ���ڻ�����ϵ�е�λ��
            %'offset',o;Ĭ����ʼ�ؽڽǣ����ڵ���ʹ��ģ����ʵ��Ƕ�һ��(defult 0)
            %'flip',logic;Ϊtrueʱ�ؽ������෴��falseһ��(default false)
            %'qlim',[a,b];�ؽڽǶ����ƣ�1X2,default[-pi,pi]��
            %'qdlim',qd;�ؽڽ��ٶ����ƣ�default[]��
            %'qddlim',qd;�ؽڽǼ��ٶ����ƣ�default[]��
            %'m'��m;����������default 0��
            %'rc',rc;������������ڻ�����ϵ��λ��������default []��
            %'I',i;��������������ĵĹ��Ծ���(default [])
            % 'Jm',J  ;���ת������(default 0)
            % 'B',B  ;�ؽ�ճ��Ħ��ϵ����default 0��
            %'torlim',tor %�ؽ���������
            % 'Tc',T  ;�������Ħ��ϵ�� (1x1 or 2x1), (default [0 0])
            % 'G',G ��������ٱ�(default 101)
            if nargin==0
                %����Ĭ������
                
                %�˶�ѧ����
                jl.w=[0 0 1]';
                jl.r=[0 0 0]';
                jl.offset=0;
                jl.flip=false;
                jl.name='';
                jl.qlim=[-pi pi];
                jl.qdlim=pi/2;%rad/s
                jl.qddlim=0.5;%rad/s^2
                jl.jointtype='R';
                %����ѧ����
                jl.m=0;
                jl.rc=[];
                jl.I=[];
                %�������
                jl.Jm=0;
                jl.B=0;
                jl.torlim=0;
                jl.Tc=0;
                jl.G=101;
            elseif nargin==1 && isa(varargin{1},'JoLink')  %��Ҫ��ȡ��Ԫ�������е�����
                %���ƶ���
                this = varargin{1};
                for j=1:length(this)
                    jl(j) = JoLink();
                    %��������
                    p = properties(this(j));
                    for i = 1:length(p)
                        jl(j).(p{i}) = this(j).(p{i});
                    end
                end
            else
                %���ݸ���������������
                %�趨��ز���
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
                    %�ж���������
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
        end%���캯��
        
        function tor=friction(jl,qd)
            %�ؽ�Ħ�������
            %�ο�Link��ͬ��������д
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
        end % Ħ����
        
        function logic=islimt(jl,q,qd,qdd,tor)
            %�����жϹؽڽǶȣ��ؽڽ��ٶȣ��Ǽ��ٶȺ����Ƿ񳬳�����
        end
    end %method
end %��Ľ���