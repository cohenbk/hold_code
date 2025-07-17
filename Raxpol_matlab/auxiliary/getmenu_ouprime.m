    while(options==process_option)
        range_km=sqrt((sum(lims(1:2))/2)^2+(sum(lims(3:4))/2)^2);
        height_km=sqrt(range_km^2+(ae)^2+2*range_km*ae*sin(ele*pi/180))-ae+0.020;
        figure(100)
        set(100,'Position',pos6)
        h11 = uicontrol('Style', 'pushbutton', 'String', 'Go up 1 tilt','Position', [2 62 100 50], 'Callback', 'options=1; uiresume');
        h21 = uicontrol('Style', 'pushbutton', 'String', 'Go down 1 tilt','Position', [2 2 100 50], 'Callback', 'options=2; uiresume');
        h61 = uicontrol('Style', 'pushbutton', 'String', 'Change Window','Position', [102 2 100 50], 'Callback', 'options=6; uiresume');
        h71 = uicontrol('Style', 'pushbutton', 'String', 'Data Cursor','Position', [102 62 100 50], 'Callback', 'options=7; uiresume');
        h81 = uicontrol('Style', 'pushbutton', 'String', 'End','Position', [302 2 100 50], 'Callback', 'options=9; uiresume');
        h1000=uicontrol(100,'Style','Text','String',['Height: ' num2str(roundn(height_km,-3)) ' km'],'Position',[15 115 130 15]);
        h1001=uicontrol(100,'Style','Text','String',['Range: ' num2str(roundn(sqrt(nanmean(lims(1:2))^2+nanmean(lims(3:4))^2),-1)) ' km'],'Position',[145 115 100 15]);
        figure(100)
        uiwait;
        close(100)
        
        if(options==1)
            ele_choice=ele_choice+1;
            if(ele_choice<max(size(ele_save{scan_num})))
                ele_number=ele_save{scan_num}(ele_choice);
            else
                scan_num=scan_num+1;
                ele_number=ele_save{scan_num}(1);
                ele_choice=1;
            end
            [row]=find(scan_number==scan_num);
            time_choice=file_time_save{scan_num}(1:6);
            filelist=files(row).name;
        elseif(options==2)
            ele_choice=ele_choice-1;
            if(ele_choice>=1)
                ele_number=ele_save{scan_num}(ele_choice);
            else
                scan_num=scan_num-1;
                ele_number=ele_save{scan_num}(max(size(ele_save{scan_num})));
                ele_choice=max(size(ele_save{scan_num}));
            end
            [row]=find(scan_number==scan_num);
            time_choice=file_time_save{scan_num}(1:6);
            filelist=files(row).name;
       elseif(options==6)
%             figure(gcf)
%             fignum=inputdlg('Figure to match limits: 1. Z, 2. Vr, 3. rhv 4. zdr'); fignum=str2double(fignum);
%             if(fignum==1)
%                 fignum=h1; fignum1=115;
%             elseif(fignum==2)
%                 fignum=h2; fignum1=3;
%             elseif(fignum==3)
%                 fignum=h3; fignum1=116;
%             else
%                 fignum=h4; fignum1=4;
%             end
%             figure(fignum1)
            lims(1:2)=get(h1,'xlim');
            lims(3:4)=get(h1,'ylim');
            figure(101)
            set(101,'Position',pos6+[0 0 50 0])
            h11 = uicontrol('Style', 'pushbutton', 'String', 'Zoom Out','Position', [2 62 100 50], 'Callback', 'options2=1; uiresume');
            h21 = uicontrol('Style', 'pushbutton', 'String', 'Zoom In','Position', [2 2 100 50], 'Callback', 'options2=2; uiresume');
            h31 = uicontrol('Style', 'pushbutton', 'String', 'Y+','Position', [102 62 100 50], 'Callback', 'options2=3; uiresume');
            h41 = uicontrol('Style', 'pushbutton', 'String', 'Y-','Position', [102 2 100 50], 'Callback', 'options2=4; uiresume');
            h51 = uicontrol('Style', 'pushbutton', 'String', 'X+','Position', [202 62 100 50], 'Callback', 'options2=5; uiresume');
            h61 = uicontrol('Style', 'pushbutton', 'String', 'X-','Position', [202 2 100 50], 'Callback', 'options2=6; uiresume');
            h71 = uicontrol('Style', 'pushbutton', 'String', 'Chg Zoom Factor','Position', [302 62 100 50], 'Callback', 'options2=7; uiresume');
            h81 = uicontrol('Style', 'pushbutton', 'String', 'Match Windows','Position', [302 2 100 50], 'Callback', 'options2=8; uiresume');
            h91 = uicontrol('Style', 'pushbutton', 'String', 'End','Position', [402 2 50 50], 'Callback', 'options2=9; uiresume');
            uiwait;
            while(options2~=9)
                if(options2==1)
                    lims(1)=lims(1)-zoom_factor_x; lims(2)=lims(2)+zoom_factor_x; lims(3)=lims(3)-zoom_factor_y; lims(4)=lims(4)+zoom_factor_y;
                elseif(options2==2)
                    lims(1)=lims(1)+zoom_factor_x; lims(2)=lims(2)-zoom_factor_x; lims(3)=lims(3)+zoom_factor_y; lims(4)=lims(4)-zoom_factor_y;
                elseif(options2==3)
                    lims(3:4)=lims(3:4)+zoom_factor_y;
                elseif(options2==4)
                    lims(3:4)=lims(3:4)-zoom_factor_y;
                elseif(options2==5)
                    lims(1:2)=lims(1:2)+zoom_factor_x;
                elseif(options2==6)
                    lims(1:2)=lims(1:2)-zoom_factor_x;
                elseif(options2==7)
                    zoom_factor_x=input(['Input new zoom factor for x (current value: ' num2str(zoom_factor_x) '): ']);
                    zoom_factor_y=input(['Input new zoom factor for y (current value: ' num2str(zoom_factor_y) '): ']);
                elseif(options2==8)
                    pause; fprintf('Change figure zoom now');
                    if(exist('pos3','var') && big_monitor==0)
                        figure(gcf)
                        fignum=inputdlg('Figure to match limits: 1. Z, 2. Vr, 3. rhv 4. zdr'); fignum=str2double(fignum);
                         figure(gcf)
                        fignum=str2double(inputdlg('Figure 1 or 2?'));
                        if(fignum==1)
                            fignum=h1;
                            fignum1=1;
                        elseif(fignum==2)
                            fignum=h2;
                            fignum1=2;
                        elseif(fignum==3)
                            fignum=h3;
                            fignum1=3;
                        elseif(fignum==4)
                            fignum=h4;
                            fignum1=4;
                        end 
                    else
                        figure(gcf)
                        fignum=str2double(inputdlg('Figure 1 or 2?'));
                        if(fignum==1)
                            fignum=h1;
                            fignum1=1;
                        else
                            fignum=h2;
                            fignum1=2;
                        end 
                    end
                    figure(fignum1)
                    lims(1:2)=get(fignum,'xlim');
                    lims(3:4)=get(fignum,'ylim');
                    zoom_factor_x=(lims(2)-lims(1))/2*0.2;
                    zoom_factor_y=(lims(4)-lims(3))/2*0.2;
                end
                    figure(1)
                    axis(lims)
                    figure(2)
                    axis(lims)
                    if(big_monitor)
                        figure(3)
                        axis(lims)
                        figure(4)
                        axis(lims)
                    end
%                 end
                figure(101)
                set(101,'Position',pos6+[0 0 50 0])
                h11 = uicontrol('Style', 'pushbutton', 'String', 'Zoom Out','Position', [2 62 100 50], 'Callback', 'options2=1; uiresume');
                h21 = uicontrol('Style', 'pushbutton', 'String', 'Zoom In','Position', [2 2 100 50], 'Callback', 'options2=2; uiresume');
                h31 = uicontrol('Style', 'pushbutton', 'String', 'Y+','Position', [102 62 100 50], 'Callback', 'options2=3; uiresume');
                h41 = uicontrol('Style', 'pushbutton', 'String', 'Y-','Position', [102 2 100 50], 'Callback', 'options2=4; uiresume');
                h51 = uicontrol('Style', 'pushbutton', 'String', 'X+','Position', [202 62 100 50], 'Callback', 'options2=5; uiresume');
                h61 = uicontrol('Style', 'pushbutton', 'String', 'X-','Position', [202 2 100 50], 'Callback', 'options2=6; uiresume');
                h71 = uicontrol('Style', 'pushbutton', 'String', 'Chg Zoom Factor','Position', [302 62 100 50], 'Callback', 'options2=7; uiresume');
                h81 = uicontrol('Style', 'pushbutton', 'String', 'Match Windows','Position', [302 2 100 50], 'Callback', 'options2=8; uiresume');
                h91 = uicontrol('Style', 'pushbutton', 'String', 'End','Position', [402 2 50 50], 'Callback', 'options2=9; uiresume');
                figure(101)
                uiwait;
                zoom_factor_x=(lims(2)-lims(1))/2*0.2;
                zoom_factor_y=(lims(4)-lims(3))/2*0.2;
%                 h1000=uicontrol(100,'Style','Text','String',['Height: ' num2str(roundn(height_km,-3)) ' km'],'Position',[135 115 130 15]);
%                 h1001=uicontrol(100,'Style','Text','String',['Range: ' num2str(roundn(sqrt(nanmean(lims(1:2))^2+nanmean(lims(3:4))^2),-1)) ' km'],'Position',[265 115 100 15]);
            end
            options=process_option;
        elseif(options==7)
            getvalue_flag=1; pdx=1; qdx=1;
            while(getvalue_flag)
                if(pdx==1)
                    pause; pdx=pdx+1;
                end
                if(az(5)-az(4)>0)
                    switch_factor=1;
                else
                    switch_factor=-1;
                end
                [x1,y1]=ginput(1);
                ran1=sqrt(x1^2+y1^2)*1000; gat1=ceil((ran1-data.detr)/data.detr); azi1=atan2(y1,x1)*180/pi;
                azi1=450-azi1-0.25*switch_factor;
                height_km=sqrt(ran1^2+(ae)^2+2*ran1*ae*sin(ele*pi/180))-ae+0.020;
                if(azi1>360)
                    azi1=azi1-360;
                end
                [c1,azi_id]=min(abs(azi1-az));
                fprintf(['Z: ' num2str(data.Z_H(azi_id,gat1)) '\n']);
                if(vr_override)
                    fprintf(['Vr: ' num2str(data.vr_h(azi_id,gat1)) '\n']);
                else
                    if(exist('vel','var'))
                        fprintf(['Vr: ' num2str(vel(gat1,azi_id)) '\n']);
                    else
                        fprintf(['Vr: ' num2str(data.vr_h(azi_id,gat1)) '\n']);
                    end
                end
                fprintf(['rhv: ' num2str(data.rho_hv(azi_id,gat1)) '\n']);
                fprintf(['Zdr: ' num2str(data.Z_DR(azi_id,gat1)) '\n']);            
%                 fprintf(['SW: ' num2str(data.SW(azi_id,gat1)) '\n']); 
                fprintf(['Height above ground: ' num2str(roundn(height_km,-3)) ' km' '\n']);
                getvalue_flag=input('Get more values? (1/0)');
            end
            options=8;
        elseif(options==9)
            is_running=false;
        end
        %Adjust limits for storm motion
        if(options==1 || options==2 || options==3 || options==4)
%             if(exist('pos3','var') && big_monitor==1)
%                 figure(115)
%                 axis(lims)
%                 figure(3)
%                 axis(lims)
%                 figure(116)
%                 axis(lims)
%                 figure(4)
%                 axis(lims)
%             else
                figure(1)
                axis(lims)
                figure(2)
                axis(lims)
                if(big_monitor)
                    figure(3)
                    axis(lims)
                    figure(4)
                    axis(lims)
                end
%             end
        end
        height_km=sqrt(range_km^2+(ae)^2+2*range_km*ae*sin(ele*pi/180))-ae+0.020;
        fprintf(['Height above ground: ' num2str(roundn(height_km,-3)) ' km' '\n']);  
    end
    