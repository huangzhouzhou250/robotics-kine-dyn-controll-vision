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


function input_button_press(h,dummy)
global robot;
global fig_1;
[filename, pathname] = uigetfile({'*.xml'},'File Selector');
robot=xml2robot3d_21(filename);
% loaddata;
intial;
end

%数据载入函数
function loaddata
global robot
n=robot.n;
setappdata(0,'F0',robot.faces{1});
setappdata(0,'P0',robot.points{1});
for i=1:n
    eval(['F',num2str(i),'=','robot.faces{i+1}',';']);
    eval(['P',num2str(i),'=','robot.points{i+1}',';']);
    setappdata(0,char(['F',num2str(i)]),robot.faces{i+1});   %此处由于第二项必须为char型数据
    setappdata(0,char(['P',num2str(i)]),robot.points{i+1}); 
end
end

%初始化函数，用以初始化界面，机器人零位显示
function intial
global robot;
global joint;
global descar;
n=robot.n;
robot.plot([ 0 0 0 0 0 0])
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
        'String','condition_number',...
        'Position',[40 150-30 20 18]);
   condition_number_edit=uicontrol(performance,'style','edit',...
        'String',[0 0 0 0 0 0],...
        'callback','@condition_number_edit_button_press',...
        'Position',[90 150-30 40 18],...
        'max',2);
   condition_number_plot=uicontrol(performance,'String','显示','callback',@condition_number_button_press,...
        'Position',[150 150-30 40 18]);
    
     %绘制零位位置
    robotplot([0 0 0 0 0 0]);

     %图形绘制的另一种方式
%     F00=getappdata(0,'F0');
%     P00=getappdata(0,'P0');
%     L0=patch('Faces',F00,'Vertices',P00);
%     set(L0,'FaceColor',[1 0 0],'EdgeColor','none');
%     H(1)=L0;
%     for i=1:n
%         temp1=getappdata(0,char(['F',num2str(i)]));
%         eval(['F',num2str(i),'0','=','temp1',';']);
%         temp2=getappdata(0,char(['P',num2str(i)]));
%         eval(['P',num2str(i),'0','=','temp2',';']);
%         temp3=patch('Faces',eval(['F',num2str(i),'0']),'Vertices',eval(['P',num2str(i),'0']));
%         eval(['L',num2str(i),'0=','temp3;']);
%         H(i+1)=eval(['L',num2str(i),'0']);
%     end
%     Tr = plot3(0,0,0,'b.');
%     H=[H,Tr];
%     setappdata(0,'patch_h',H);
%     setappdata(0,'ThetaOld',zeros(1,n));

end

%绘图函数，绘制指定位姿的下的机器人三维模型
function robotplot(theta)
%输入theta为机器人的姿态
%输入后会修改机器人各个参数
global robot;
global fig_1;
n=robot.n;
%绘制图形
robot.plot(theta);

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
global tdr_edit tdx_edit tdy_edit tdp_edit  tdyo_edit tdz_edit;
set(tdx_edit,'string',num2str(p(1)));
set(tdy_edit,'string',num2str(p(2)));
set(tdz_edit,'string',num2str(p(3)));
set(tdr_edit,'string',num2str(rpy(1)));
set(tdp_edit,'string',num2str(rpy(2)));
set(tdyo_edit,'string',num2str(rpy(3)));


end
%按钮函数
function tdx_edit_button_press(h,dummy)

end
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

