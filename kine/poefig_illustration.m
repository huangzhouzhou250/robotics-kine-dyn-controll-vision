%% poefig��غ�����¼��˵��
%condition_number_button_press ��������ʾ��ť  ������ͼ������ʾĳ�����켣��������
%condition_number_edit_button_press �������༭�ؼ�
%manipulability_button_press �ɲ�������ʾ��ť ������ͼ������ʾĳ�����켣��������
%manipulability_edit_button_press �ɲ�������ʾ��ť
%screwi_edit_button_press ��i����������ʾ��ť
%inputjoint_button_press ����Ƕȹ켣������ť
%inputdescar_button_press ����ĩ��ִ�����켣������ť
%tci_edit_button_press  �ο��Ƕȱ༭�ؼ�
%tdx_edit_button_press,tdy_edit_button_press,tdz_edit_button_pressĩ��ִ����λ�ã�x,y,z���༭�ؼ�
%tdr_edit_button_press��tdp_edit_button_press��tdyo_edit_button_pressĩ��ִ�����Ƕȣ�rpy���༭�ؼ�
%ti_slider_button_press��ti_edit_button_press �ؽ�i�ǶȻ����ؼ��ͱ༭�ؼ�

%% ���ݴ洢˵��
% �洢�켣λ��
% setappdata(0,'xtrail',0);
% setappdata(0,'ytrail',0);
% setappdata(0,'ztrail',0)

% setappdata(0,'Tr',Tr); %ĩ��λ�û�ͼ

% setappdata(0,'fig0',fig); %��ͼָ��������

% setappdata(0,'ThetaOld',theta); %��һ��״̬�ĹؽڽǶ�

%% ��Ҫ����Ĺ���
%1.���ݻ����˵Ĵ�С��ͼ��������

%2.������ʾ���˶�ѧ����δ��ʾ

%3.״̬�ռ䵼�������ģ��


%����˵��
     %ͼ�λ��Ƶ���һ�ַ�ʽ
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