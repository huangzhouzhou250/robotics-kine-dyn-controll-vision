%基于指数积公式的机器人运动学仿真平台

%% 绘制基础框架
%绘制框架fig
set(0,'Units','pixels')
dim = get(0,'ScreenSize');
global fig_1;
fig_1 = figure('doublebuffer','on','Position',[0,35,dim(3)-150,dim(4)-100],...
    'MenuBar','none','Name','POE Robot Drawing',...
    'NumberTitle','off');
hold on;

light    %灯光                           
daspect([1 1 1])                   
view(135,25)   %视角
xlabel('X'),ylabel('Y'),zlabel('Z');
title('指数积机器人仿真','FontSize',14);
axis([-800 800 -800 800 -100 1100]);
%绘制图形边框 
plot3([-800,800],[-800,-800],[-100,-100],'k')
plot3([-800,-800],[-800,800],[-100,-100],'k')
plot3([-800,-800],[-800,-800],[-100,1100],'k')
plot3([-800,800],[-800,-800],[1100,1100],'k')
plot3([-800,-800],[-800,800],[1100,1100],'k')
plot3([-800,-800],[800,800],[-100,1100],'k')
Tr=plot3(0,0,0,'b.');
setappdata(0,'xtrail',0);
setappdata(0,'ytrail',0);
setappdata(0,'ztrail',0);
setappdata(0,'Tr',Tr);

fig=fig_1.Children(1);
setappdata(0,'fig0',fig);

global robot;

%绘制主要功能区 
global joint;
joint = uipanel(fig_1,...
    'Position',[0.02 0.05  0.23 0.4],...
    'Title','关节角度','FontSize',11);

global descar;
descar=uipanel(fig_1,...
    'Position',[0.02 0.5 0.23 0.2],...
    'Title','末端位姿','FontSize',11);

func=uipanel(fig_1,...
    'Position',[0.02 0.9 0.3 0.1],...
     'FontSize',11,'BorderType','none');
 
global tracj_input;
tracj_input=uipanel(fig_1,...
    'Position',[0.02 0.72 0.23 0.13],...
    'Title','轨迹导入','FontSize',11);

global screw_display;
screw_display=uipanel(fig_1,...
    'Position',[0.8 0.1 0.2 0.4],...
    'Title','实时旋量','FontSize',11);

global performance;
performance=uipanel(fig_1,...
    'Position',[0.8 0.55 0.2 0.3],...
    'Title','运动学性能','FontSize',11);

input = uicontrol(func,'String','导入','callback',@input_button_press,...
    'Position',[0 20 30 30]);

trail_delete = uicontrol(func,'String','删除轨迹','callback',@trail_delete_button_press,...
    'Position',[60 20 50 30]);

home = uicontrol(func,'String','初始位置','callback',@home_button_press,...
    'Position',[120 20 50 30]);

% 导入机器人的xml文件，并初始化相应界面
function input_button_press(h,dummy)
global robot;

%读取机器人的xml文件，并生成对应的机器人模型
[filename, pathname] = uigetfile({'*.xml'},'File Selector'); %读取机器人xml文件
robot=xml2robot3d_21(filename);%生成机器人函数

% 初始化界面并记录关节零位位置
intial;

%获取并记录初始位置
q0=robot.offset;
setappdata(0,'ThetaOld',q0);
% robotanimation([1 1 1 1 1 1],50,'y'); %用于检查动画函数是否正确
end

%初始化函数，用以初始化界面，机器人零位显示
function intial
    global robot;
    global joint;
    global descar;
    n=robot.n;
    %绘制关节角度控制按钮
    LD = 105; % Left, used to set the GUI.
    HT = 18;  % Height
    BT = 210; % Bottom
    qlim=robot.qlim;
    for i=1:n
        temp_sl=uicontrol(joint,'style','slider',...
            'Max',qlim(i,2),'Min',qlim(i,1),'Value',0,...
            'SliderStep',[0.05 0.2],...
            'callback',eval(['@t',num2str(i),'_slider_button_press']),...
            'Position',[LD BT-30*(i-1) 120 HT]);
        eval(['global t',num2str(i),'_slider;']);
        eval(['t',num2str(i),'_slider=temp_sl;']);
        temp_min=uicontrol(joint,'style','text',...
            'String',string(qlim(i,1)),...
            'Position',[LD-30 BT+1-30*(i-1) 25 HT-4]);
        eval(['t',num2str(i),'_min=temp_min;']);
        temp_max= uicontrol(joint,'style','text',...
            'String',string(qlim(i,2)),...
            'Position',[LD+125 BT+1-30*(i-1) 30 HT-4]);
        eval(['t',num2str(i),'_max=temp_max']);
        temp_text=uibutton(joint,'style','text',... 
            'String',['\theta_',num2str(i)],...               
            'Position',[LD-100 BT-30*(i-1) 20 HT]); 
        eval(['t',num2str(i),'_text=temp_text;']);
        temp_edit=uicontrol(joint,'style','edit',...
            'String',0,...
            'callback',eval(['@t',num2str(i),'_edit_button_press']),...
            'Position',[LD-75 BT-30*(i-1) 30 HT]); % L, B, W, H
        eval(['global t',num2str(i),'_edit;']);
        eval(['t',num2str(i),'_edit=temp_edit;']);
    end
    
    %绘制笛卡尔空间(逆运动学)控制按钮
    %末端执行器角度绘制
    ikine = uicontrol(descar,'String','逆运动学','callback',@ikine_button_press,...
        'Position',[115 10 50 20]);
    jiaodu_text = uibutton(descar,'style','text',...
        'String','角度(deg)','FontSize',10,...
        'Position',[7 90 50 18]); 
    tdr_text = uibutton(descar,'style','text',...
        'String','r',...
        'Position',[10 70 30 18]);
    global tdr_edit tdp_edit tdyo_edit tdx_edit tdy_edit tdz_edit;
    tdr_edit = uicontrol(descar,'style','edit',...
        'String',0,...
        'callback',@tdr_edit_button_press,...
        'Position',[30 70 30 18]);
    tdp_text = uibutton(descar,'style','text',...
        'String','p',...
        'Position',[10 45 30 18]);
    tdp_edit = uicontrol(descar,'style','edit',...
        'String',0,...
        'callback',@tdp_edit_button_press,...
        'Position',[30 45 30 18]);
    tdyo_text = uibutton(descar,'style','text',...
        'String','y',...
        'Position',[10 20 30 18]);
    tdyo_edit = uicontrol(descar,'style','edit',...
        'String',0,...
        'callback',@tdyo_edit_button_press,...
        'Position',[30 20 30 18]);
    %末端执行器位置绘制
    weizhi_text = uibutton(descar,'style','text',...
        'String','位置(mm)','FontSize',10,...
        'Position',[65 90 50 18]);
    tdx_text = uibutton(descar,'style','text',...
        'String','x',...
        'Position',[65 70 30 18]);
    tdx_edit = uicontrol(descar,'style','edit',...
        'String',0,...
        'callback',@tdx_edit_button_press,...
        'Position',[85 70 30 18]);
    tdy_text = uibutton(descar,'style','text',...
        'String','y',...
        'Position',[65 45 30 18]);
    tdy_edit = uicontrol(descar,'style','edit',...
        'String',0,...
        'callback',@tdy_edit_button_press,...
        'Position',[85 45 30 18]);
    tdz_text = uibutton(descar,'style','text',...
        'String','z',...
        'Position',[65 20 30 18]);
    tdz_edit = uicontrol(descar,'style','edit',...
        'String',0,...
        'callback',@tdz_edit_button_press,...
        'Position',[85 20 30 18]);

    %参考角度控件绘制
    cankaojiaodu_text = uibutton(descar,'style','text',...
    'String','参考角度（deg）','FontSize',10,...
    'Position',[195 90 50 18]); 
    for i=1:n
        if i<=3
            temp_text=uibutton(descar,'style','text',...
                'String',['\theta_',num2str(i)],...
                'Position',[165 70-25*(i-1) 30 18]);
            eval(['tc',num2str(i),'_text=temp_text;']);
            temp_edit=uicontrol(descar,'style','edit',...
                'String',0,...
                'callback',eval(['@tc',num2str(i),'_edit_button_press']),...
                'Position',[185 70-25*(i-1) 30 18]);
            eval(['tc',num2str(i),'_edit=temp_edit;']); 
        else
            temp_text=uibutton(descar,'style','text',...
                'String',['\theta_',num2str(i)],...
                'Position',[210 70-25*(i-4) 30 18]);
            eval(['tc',num2str(i),'_text=temp_text;']);
            temp_edit=uicontrol(descar,'style','edit',...
                'String',0,...
                'callback',eval(['@tc',num2str(i),'_edit_button_press']),...
                'Position',[235 70-25*(i-4) 30 18]);
            eval(['tc',num2str(i),'_edit=temp_edit;']); 
        end
    end
    
    %轨迹输入控件绘制
    global tracj_input;
    input_joint = uicontrol(tracj_input,'String','角度轨迹','callback',@inputjoint_button_press,...
        'Position',[50 20 50 30]);
    input_descar = uicontrol(tracj_input,'String','末端轨迹','callback',@inputdescar_button_press,...
        'Position',[150 20 50 30]);
    
    %旋量输出控件绘制
    global screw_display;
    for i=1:n
        temp_text= uibutton(screw_display,'style','text',...
            'String',['\xi_',num2str(i)],...
            'Position',[10 220-25*(i-1) 20 18]);
        eval(['twist',num2str(i),'_text=temp_text;']);
        temp_edit=uicontrol(screw_display,'style','edit',...
            'String',[0 0 0 0 0 0],...
            'callback',eval(['@screw',num2str(i),'_edit_button_press']),...
            'Position',[30 220-25*(i-1) 180 18],...
            'max',2);
        eval(['twist',num2str(i),'_edit=temp_edit;']);
    end
    
    %运动学性能控件绘制
    global performance;
    manipulability_text=uibutton(performance,'style','text',...
        'String','manipulability',...
        'Position',[40 150 20 18]);
    manipulability_edit=uicontrol(performance,'style','edit',...
        'String',[0 0 0 0 0 0],...
        'callback','@manipulability_edit_button_press',...
        'Position',[90 150 40 18],...
        'max',2);
    manipulability_plot=uicontrol(performance,'String','显示','callback',@manipulability_button_press,...
        'Position',[150 150 40 18]);

    condition_number_text=uibutton(performance,'style','text',...
        'String','条件数',...
        'Position',[40 150-30 20 18]);
   condition_number_edit=uicontrol(performance,'style','edit',...
        'String',[0 0 0 0 0 0],...
        'callback','@condition_number_edit_button_press',...
        'Position',[90 150-30 40 18],...
        'max',2);
   condition_number_plot=uicontrol(performance,'String','显示','callback',@condition_number_button_press,...
        'Position',[150 150-30 40 18]);
    
     %绘制零位位置
    robotplot([robot.offset],'n');
end

%绘图函数，绘制指定位姿的下的机器人三维模型
function robotplot(theta,trail)
%输入theta为机器人的姿态
%输入后会修改机器人各个参数，并绘制机器人的三维图形
if nargin==1
    trail='n';
elseif nargin==2
    if ~strcmp(trail,'n') && ~strcmp(trail,'y')
        error('轨迹类型输入有误')
    end
else
    error('输入有误');
end

global robot;
n=robot.n;
%绘制图形
fig=getappdata(0,'fig0');
axes(fig);
child=fig.Children;
child_length=length(child);
if child_length>8
    delete(child(1:child_length-8));
end
h=robot.plot(theta);

%修改角度框参数
q0=theta*180/pi;
for i=1:n
    eval(['global t',num2str(i),'_edit;']);
    eval(['global t',num2str(i),'_slider;']);
    set(eval(['t',num2str(i),'_edit']),'string',num2str(q0(i)));
    set(eval(['t',num2str(i),'_slider']),'value',q0(i));
end

%修改末端执行器和位置
TE=robot.fkinep(theta);
[R,p]=tr2rt(TE);
rpy=tr2rpy(R);
rpy=rpy*180/pi;
global tdr_edit tdx_edit tdy_edit tdp_edit  tdyo_edit tdz_edit;
set(tdx_edit,'string',num2str(p(1)));
set(tdy_edit,'string',num2str(p(2)));
set(tdz_edit,'string',num2str(p(3)));
set(tdr_edit,'string',num2str(rpy(1)));
set(tdp_edit,'string',num2str(rpy(2)));
set(tdyo_edit,'string',num2str(rpy(3)));

if trail=='y'
    x_trail = getappdata(0,'xtrail');
    y_trail = getappdata(0,'ytrail');
    z_trail = getappdata(0,'ztrail');
    xdata = [x_trail p(1)];
    ydata = [y_trail p(2)];
    zdata = [z_trail p(3)];
    setappdata(0,'xtrail',xdata); 
    setappdata(0,'ytrail',ydata); 
    setappdata(0,'ztrail',zdata);
    Tr=getappdata(0,'Tr');
    set(Tr,'xdata',xdata,'ydata',ydata,'zdata',zdata);
end
end

%绘图函数，绘制机器人的动画并改变相应控件值
function robotanimation(theta,n,trail,method)
%theta为机器人的目标关节角度值
%n为插入的步数
%trail用于确定是否显示轨迹
%trail='n' 不显示轨迹
%trail='y' 显示轨迹
%method用于确定机器人起始点与目标点之间的插值方式
%1 使用五次插值函数(traj_5)
%2 使用三次插值函数(tracj_3)
%3 使用梯形插值函数(tracj_t)

global robot;

%判断输入变量个数
if nargin==1
    n=30;
    trail='n';
    method=1;
elseif nargin==2
    trail='n';
    method=1;
elseif nargin==3
    method=1;
    if ~strcmp(trail,'n')&&~strcmp(trail,'y')
        error('trail类型输入有误')
    end
elseif nargin==4
    if method~=1||method~=2||method~=3
        error('method类型输入有误')
    end
end

%判断输入角度长度是否正确
n_robot=robot.n;
if length(theta)~=n_robot
    error('角度输入长度有误')
end

%读取之前角度
theta0=getappdata(0,'ThetaOld');

if method==1
    q=traj_5(theta0,theta,n);
elseif method==2
    q=tracj_3(theta0,theta,n);
elseif method==3
    qdd=robot.qddlim;
    q=tracj_t(theta0,theta,5,'A',qdd,n);
end
n_plot=size(q,1);

for i=2:1:n_plot
    q_temp=q(i,:);
    robotplot(q_temp,trail);
    drawnow;
end
setappdata(0,'ThetaOld',theta);
end


%关节1修改函数
function t1_slider_button_press(h,dummy)
global t1_edit;
slider_value = round(get(h,'Value'));
set(t1_edit,'string',slider_value);
q0=getappdata(0,'ThetaOld');
q0(1)=slider_value*pi/180;
robotanimation(q0,20,'n');
end

function t1_edit_button_press(h,dummy)
global robot t1_edit t1_slider;
qlim=robot.qlim;
user_entry = check_edit(h,qlim(1,1),qlim(1,2),0,t1_edit);
set(t1_slider,'Value',user_entry);
q0=getappdata(0,'ThetaOld');
q0(1)=user_entry*pi/180;
robotanimation(q0,20,'n')
end

%关节2修改函数
function t2_slider_button_press(h,dummy)
global t2_edit;
slider_value = round(get(h,'Value'));
set(t2_edit,'string',slider_value);
q0=getappdata(0,'ThetaOld');
q0(2)=slider_value*pi/180;
robotanimation(q0,20,'n');
end

function t2_edit_button_press(h,dummy)
global robot t2_edit t2_slider;
qlim=robot.qlim;
user_entry = check_edit(h,qlim(2,1),qlim(2,2),0,t2_edit);
set(t2_slider,'Value',user_entry);
q0=getappdata(0,'ThetaOld');
q0(2)=user_entry*pi/180;
robotanimation(q0,20,'n')
end

%关节3修改函数
function t3_slider_button_press(h,dummy)
global t3_edit;
slider_value = round(get(h,'Value'));
set(t3_edit,'string',slider_value);
q0=getappdata(0,'ThetaOld');
q0(3)=slider_value*pi/180;
robotanimation(q0,20,'n');
end

function t3_edit_button_press(h,dummy)
global robot t3_edit t3_slider;
qlim=robot.qlim;
user_entry = check_edit(h,qlim(3,1),qlim(3,2),0,t3_edit);
set(t3_slider,'Value',user_entry);
q0=getappdata(0,'ThetaOld');
q0(3)=user_entry*pi/180;
robotanimation(q0,20,'n')
end

%关节4修改函数
function t4_slider_button_press(h,dummy)
global t4_edit;
slider_value = round(get(h,'Value'));
set(t4_edit,'string',slider_value);
q0=getappdata(0,'ThetaOld');
q0(4)=slider_value*pi/180;
robotanimation(q0,20,'n');
end

function t4_edit_button_press(h,dummy)
global robot t4_edit t4_slider;
qlim=robot.qlim;
user_entry = check_edit(h,qlim(4,1),qlim(4,2),0,t4_edit);
set(t4_slider,'Value',user_entry);
q0=getappdata(0,'ThetaOld');
q0(4)=user_entry*pi/180;
robotanimation(q0,20,'n')
end

%关节5修改函数
function t5_slider_button_press(h,dummy)
global t5_edit;
slider_value = round(get(h,'Value'));
set(t5_edit,'string',slider_value);
q0=getappdata(0,'ThetaOld');
q0(5)=slider_value*pi/180;
robotanimation(q0,20,'n');
end

function t5_edit_button_press(h,dummy)
global robot t5_edit t5_slider;
qlim=robot.qlim;
user_entry = check_edit(h,qlim(5,1),qlim(5,2),0,t1_edit);
set(t5_slider,'Value',user_entry);
q0=getappdata(0,'ThetaOld');
q0(5)=user_entry*pi/180;
robotanimation(q0,20,'n')
end

%关节6修改函数
function t6_slider_button_press(h,dummy)
global t6_edit;
slider_value = round(get(h,'Value'));
set(t6_edit,'string',slider_value);
q0=getappdata(0,'ThetaOld');
q0(6)=slider_value*pi/180;
robotanimation(q0,20,'n');
end

function t6_edit_button_press(h,dummy)
global robot t6_edit t6_slider;
qlim=robot.qlim;
user_entry = check_edit(h,qlim(6,1),qlim(6,2),0,t6_edit);
set(t6_slider,'Value',user_entry);
q0=getappdata(0,'ThetaOld');
q0(6)=user_entry*pi/180;
robotanimation(q0,20,'n')
end

%关节7修改函数
function t7_slider_button_press(h,dummy)
global t7_edit;
slider_value = round(get(h,'Value'));
set(t7_edit,'string',slider_value);
q0=getappdata(0,'ThetaOld');
q0(7)=slider_value*pi/180;
robotanimation(q0,20,'n');
end

function t7_edit_button_press(h,dummy)
global robot t7_edit t7_slider;
qlim=robot.qlim;
user_entry = check_edit(h,qlim(7,1),qlim(7,2),0,t7_edit);
set(t7_slider,'Value',user_entry);
q0=getappdata(0,'ThetaOld');
q0(7)=user_entry*pi/180;
robotanimation(q0,20,'n')
end

%关节1修改函数
function t8_slider_button_press(h,dummy)
global t8_edit;
slider_value = round(get(h,'Value'));
set(t8_edit,'string',slider_value);
q0=getappdata(0,'ThetaOld');
q0(81)=slider_value*pi/180;
robotanimation(q0,20,'n');
end

function t8_edit_button_press(h,dummy)
global robot t8_edit t8_slider;
qlim=robot.qlim;
user_entry = check_edit(h,qlim(8,1),qlim(8,2),0,t8_edit);
set(t8_slider,'Value',user_entry);
q0=getappdata(0,'ThetaOld');
q0(8)=user_entry*pi/180;
robotanimation(q0,20,'n')
end

%修改末端位置函数 
function tdx_edit_button_press(h,dummy)

end

%轨迹删除函数
function trail_delete_button_press(h,dummy)

Tr=getappdata(0,'Tr');%获取轨迹

%将轨迹的值设置为0
setappdata(0,'xtrail',0);
setappdata(0,'ytrail',0);
setappdata(0,'ztrail',0);
%清除轨迹
set(Tr,'xdata',0,'ydata',0,'zdata',0);
end

%回到零位位置
function home_button_press(h,dummy)
global robot;
robotplot(robot.offset);
Tr=getappdata(0,'Tr');%获取轨迹

%将轨迹的值设置为0
setappdata(0,'xtrail',0);
setappdata(0,'ytrail',0);
setappdata(0,'ztrail',0);
%清除轨迹
set(Tr,'xdata',0,'ydata',0,'zdata',0);
setappdata(0,'ThetaOld',robot.offset);
end

%控件生成函数
function [hout,ax_out] = uibutton(varargin)
        %uibutton: Create pushbutton with more flexible labeling than uicontrol.
        % Usage:
        %   uibutton accepts all the same arguments as uicontrol except for the
        %   following property changes:
        %
        %     Property      Values
        %     -----------   ------------------------------------------------------
        %     Style         'pushbutton', 'togglebutton' or 'text', default =
        %                   'pushbutton'.
        %     String        Same as for text() including cell array of strings and
        %                   TeX or LaTeX interpretation.
        %     Interpreter   'tex', 'latex' or 'none', default = default for text()
        %
        % Syntax:
        %   handle = uibutton('PropertyName',PropertyValue,...)
        %   handle = uibutton(parent,'PropertyName',PropertyValue,...)
        %   [text_obj,axes_handle] = uibutton('Style','text',...
        %       'PropertyName',PropertyValue,...)
        %
        % uibutton creates a temporary axes and text object containing the text to
        % be displayed, captures the axes as an image, deletes the axes and then
        % displays the image on the uicontrol.  The handle to the uicontrol is
        % returned.  If you pass in a handle to an existing uicontol as the first
        % argument then uibutton will use that uicontrol and not create a new one.
        %
        % If the Style is set to 'text' then the axes object is not deleted and the
        % text object handle is returned (as well as the handle to the axes in a
        % second output argument).
        %
        % See also UICONTROL.

        % Version: 1.6, 20 April 2006
        % Author:  Douglas M. Schwarz
        % Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
        % Real_email = regexprep(Email,{'=','*'},{'@','.'})


        % Detect if first argument is a uicontrol handle.
        keep_handle = false;
        if nargin > 0
            h = varargin{1};
            if isscalar(h) && ishandle(h) && strcmp(get(h,'Type'),'uicontrol')
                keep_handle = true;
                varargin(1) = [];
            end
        end

        % Parse arguments looking for 'Interpreter' property.  If found, note its
        % value and then remove it from where it was found.
        interp_value = get(0,'DefaultTextInterpreter');
        arg = 1;
        remove = [];
        while arg <= length(varargin)
            v = varargin{arg};
            if isstruct(v)
                fn = fieldnames(v);
                for i = 1:length(fn)
                    if strncmpi(fn{i},'interpreter',length(fn{i}))
                        interp_value = v.(fn{i});
                        v = rmfield(v,fn{i});
                    end
                end
                varargin{arg} = v;
                arg = arg + 1;
            elseif ischar(v)
                if strncmpi(v,'interpreter',length(v))
                    interp_value = varargin{arg+1};
                    remove = [remove,arg,arg+1];
                end
                arg = arg + 2;
            elseif arg == 1 && isscalar(v) && ishandle(v) && ...
                    any(strcmp(get(h,'Type'),{'figure','uipanel'}))
                arg = arg + 1;
            else
                error('Invalid property or uicontrol parent.')
            end
        end
        varargin(remove) = [];

        % Create uicontrol, get its properties then hide it.
        if keep_handle
            set(h,varargin{:})
        else
            h = uicontrol(varargin{:});
        end
        s = get(h);
        if ~any(strcmp(s.Style,{'pushbutton','togglebutton','text'}))
            delete(h)
            error('''Style'' must be pushbutton, togglebutton or text.')
        end
        set(h,'Visible','off')

        % Create axes.
        parent = get(h,'Parent');
        ax = axes('Parent',parent,...
            'Units',s.Units,...
            'Position',s.Position,...
            'XTick',[],'YTick',[],...
            'XColor',s.BackgroundColor,...
            'YColor',s.BackgroundColor,...
            'Box','on',...
            'Color',s.BackgroundColor);
        % Adjust size of axes for best appearance.
        set(ax,'Units','pixels')
        pos = round(get(ax,'Position'));
        if strcmp(s.Style,'text')
            set(ax,'Position',pos + [0 1 -1 -1])
        else
            set(ax,'Position',pos + [4 4 -8 -8])
        end
        switch s.HorizontalAlignment
            case 'left'
                x = 0.0;
            case 'center'
                x = 0.5;
            case 'right'
                x = 1;
        end
        % Create text object.
        text_obj = text('Parent',ax,...
            'Position',[x,0.5],...
            'String',s.String,...
            'Interpreter',interp_value,...
            'HorizontalAlignment',s.HorizontalAlignment,...
            'VerticalAlignment','middle',...
            'FontName',s.FontName,...
            'FontSize',s.FontSize,...
            'FontAngle',s.FontAngle,...
            'FontWeight',s.FontWeight,...
            'Color',s.ForegroundColor);

        % If we are creating something that looks like a text uicontrol then we're
        % all done and we return the text object and axes handles rather than a
        % uicontrol handle.
        if strcmp(s.Style,'text')
            delete(h)
            if nargout
                hout = text_obj;
                ax_out = ax;
            end
            return
        end

        % Capture image of axes and then delete the axes.
        frame = getframe(ax);
        delete(ax)

        % Build RGB image, set background pixels to NaN and put it in 'CData' for
        % the uicontrol.
        if isempty(frame.colormap)
            rgb = frame.cdata;
        else
            rgb = reshape(frame.colormap(frame.cdata,:),[pos([4,3]),3]);
        end
        size_rgb = size(rgb);
        rgb = double(rgb)/255;
        back = repmat(permute(s.BackgroundColor,[1 3 2]),size_rgb(1:2));
        isback = all(rgb == back,3);
        rgb(repmat(isback,[1 1 3])) = NaN;
        set(h,'CData',rgb,'String','','Visible',s.Visible)

        % Assign output argument if necessary.
        if nargout
            hout = h;
        end
%%
end
%控件相关设置函数
function user_entry = check_edit(h,min_v,max_v,default,h_edit)
        % This function will check the value typed in the text input box
        % against min and max values, and correct errors.
        %
        % h: handle of gui
        % min_v min value to check
        % max_v max value to check
        % default is the default value if user enters non number
        % h_edit is the edit value to update.
        %
        user_entry = str2double(get(h,'string'));
        if isnan(user_entry)
            errordlg(['You must enter a numeric value, defaulting to ',num2str(default),'.'],'Bad Input','modal')
            set(h_edit,'string',default);
            user_entry = default;
        end
        %
        if user_entry < min_v
            errordlg(['Minimum limit is ',num2str(min_v),' degrees, using ',num2str(min_v),'.'],'Bad Input','modal')
            user_entry = min_v;
            set(h_edit,'string',user_entry);
        end
        if user_entry > max_v
            errordlg(['Maximum limit is ',num2str(max_v),' degrees, using ',num2str(max_v),'.'],'Bad Input','modal')
            user_entry = max_v;
            set(h_edit,'string',user_entry);
        end
    end
%
