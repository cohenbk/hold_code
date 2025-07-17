    while(options==process_option)
        range_km=sqrt((sum(lims(1:2))/2)^2+(sum(lims(3:4))/2)^2);
%         height_km=sqrt((sum(lims(1:2))/2)^2+(sum(lims(3:4))/2)^2)*sin(ele*pi/180)+0.020;
        height_km=sqrt(range_km^2+(ae)^2+2*range_km*ae*sin(ele*pi/180))-ae+0.020;
        figure(100)
        set(100,'Position',pos6)
        h11 = uicontrol('Style', 'pushbutton', 'String', 'Go up 1 tilt','Position', [2 62 100 50], 'Callback', 'options=1; uiresume');
        h21 = uicontrol('Style', 'pushbutton', 'String', 'Go down 1 tilt','Position', [2 2 100 50], 'Callback', 'options=2; uiresume');
        h31 = uicontrol('Style', 'pushbutton', 'String', 'Go forward 1 vol','Position', [102 62 100 50], 'Callback', 'options=3; uiresume');
        h41 = uicontrol('Style', 'pushbutton', 'String', 'Go back 1 vol','Position', [102 2 100 50], 'Callback', 'options=4; uiresume');
        h51 = uicontrol('Style', 'pushbutton', 'String', 'Change Time','Position', [202 62 100 50], 'Callback', 'options=5; uiresume');
        h61 = uicontrol('Style', 'pushbutton', 'String', 'Change Window','Position', [202 2 100 50], 'Callback', 'options=6; uiresume');
        h71 = uicontrol('Style', 'pushbutton', 'String', 'Data Cursor','Position', [302 62 100 50], 'Callback', 'options=7; uiresume');
        h81 = uicontrol('Style', 'pushbutton', 'String', 'End','Position', [302 2 100 50], 'Callback', 'options=9; uiresume');
        h91 = uicontrol(100,'Style', 'pushbutton', 'String', 'Chg Left Var','Position', [402 62 100 50], 'Callback', 'options=11; uiresume;');
        h101 = uicontrol(100,'Style', 'pushbutton', 'String', 'Chg Right Var','Position', [402 2 100 50], 'Callback', 'options=12; uiresume;');
        h111= uicontrol(100,'Style','pushbutton','String','Save Fig','Position',[502 2 100 50],'Callback','options=13; uiresume;');
        h112= uicontrol(100,'Style','pushbutton','String','RHI','Position',[502 62 100 50],'Callback','options=14; uiresume;');
        h1000=uicontrol(100,'Style','Text','String',['Height: ' num2str(roundn(height_km,-3)) ' km'],'Position',[135 115 130 15]);
        h1001=uicontrol(100,'Style','Text','String',['Range: ' num2str(roundn(sqrt(nanmean(lims(1:2))^2+nanmean(lims(3:4))^2),-1)) ' km'],'Position',[265 115 100 15]);
        figure(100)
        uiwait;
        close(100)
        
        if(options==11)
            if(plot_choice(1)==5)
                figure(1)
                colorbar off
            end
            new_color_flag=[1 0];
            figure(102)
            set(102,'Position',pos7)
    %         text(pos7(1),pos7(2)/2,'Left Panel Control')
            h12 = uicontrol(102,'Style', 'pushbutton', 'String', 'Z','Position',[2 2 60 40], 'Callback', 'plot_choice(1)=1; options=10; uiresume;');
            h22 = uicontrol(102,'Style', 'pushbutton', 'String', 'Vr','Position',[2 44 60 40], 'Callback', 'plot_choice(1)=2; options=10; uiresume;');
            h32 = uicontrol(102,'Style', 'pushbutton', 'String', 'rhohv','Position',[64 2 60 40], 'Callback', 'plot_choice(1)=3; options=10; uiresume;');
            h42 = uicontrol(102,'Style', 'pushbutton', 'String', 'ZDR','Position',[64 44 60 40], 'Callback', 'plot_choice(1)=4; options=10; uiresume;');
            h52 = uicontrol(102,'Style', 'pushbutton', 'String', 'Thres','Position',[126 44 60 40], 'Callback', 'plot_choice(1)=5; options=10; uiresume;');
            uiwait;
            close(102)
        elseif(options==12)
            if(plot_choice(2)==5)
                figure(2)
                colorbar off
            end
            new_color_flag=[0 1];
            figure(103)
            set(103,'Position',pos8)
    %         text(pos8(1),pos8(2)/2,'Right Panel Control')
            h13 = uicontrol(103,'Style', 'pushbutton', 'String', 'Z','Position',[2 2 60 40], 'Callback', 'plot_choice(2)=1; options=10; uiresume;');
            h14 = uicontrol(103,'Style', 'pushbutton', 'String', 'Vr','Position',[2 44 60 40], 'Callback', 'plot_choice(2)=2; options=10; uiresume;');
            h15 = uicontrol(103,'Style', 'pushbutton', 'String', 'rhohv','Position',[64 2 60 40], 'Callback', 'plot_choice(2)=3; options=10; uiresume;');
            h16 = uicontrol(103,'Style', 'pushbutton', 'String', 'ZDR','Position',[64 44 60 40], 'Callback', 'plot_choice(2)=4; options=10; uiresume;');
            h17 = uicontrol(103,'Style', 'pushbutton', 'String', 'Thres','Position',[126 44 60 40], 'Callback', 'plot_choice(2)=5; options=10; uiresume;');
            uiwait;
            close(103)
        elseif(options==13)
            if(data_cursor_mode)
                fprintf('Warning: in data cursor mode bad xx,yy values possible');
                pause;
            end
             
            title_str={'Z','vr','rhohv','ZDR','thres'};
            if(data_flags(1)==1)
                if(~exist('tornado_number','var'))
                    tornado_number=input(['Which storm? (1: El Reno, 2: Chickasha/Newcastle, 3: Goldsby): ']);
                    while(min(size(tornado_number))<1)
                        tornado_number=input(['Which storm? (1: El Reno, 2: Chickasha/Newcastle, 3: Goldsby): ']);
                    end
                end
            else
                tornado_number=1;
            end
            fig_dir=['/Users/alexlyakhov/Documents/MATLAB/figs'];
                    
            if(~exist(fig_dir))
                mkdir(fig_dir);
            end
            figure(1)
            picname=[num2str(tornado_number) '_' time_choice '_' num2str(ele_save{scan_num}(ele_choice)) '_' title_str{plot_choice(1)} '.png'];
            cd(fig_dir)
            print(picname,'-dpng');
            cd(root_dir)
            figure(2)
            picname=[num2str(tornado_number) '_' time_choice '_' num2str(ele_save{scan_num}(ele_choice)) '_' title_str{plot_choice(2)} '.png'];
            cd(fig_dir)
            print(picname,'-dpng');
            cd(root_dir)
            if(big_monitor)
                figure(3)
                picname=[num2str(tornado_number) '_' time_choice '_' num2str(ele_save{scan_num}(ele_choice)) '_' title_str{plot_choice(3)} '.png'];
                cd(fig_dir)
                print(picname,'-dpng');
                cd(root_dir)
                figure(4)
                picname=[num2str(tornado_number) '_' time_choice '_' num2str(ele_save{scan_num}(ele_choice)) '_' title_str{plot_choice(4)} '.png'];
                cd(fig_dir)
                print(picname,'-dpng');
                cd(root_dir)
                figure(5)
                picname=[num2str(tornado_number) '_' time_choice '_' num2str(ele_save{scan_num}(ele_choice)) '_' title_str{plot_choice(5)} '.png'];
                cd(fig_dir)
                print(picname,'-dpng');
                cd(root_dir)
            end
        elseif(options==14)
            plot_rhi;
        end
        
%         options=menu('Pick an option','1: Go up 1 tilt','2: Go down 1 tilt','3: Go forward 1 volume','4: Go back 1 volume','5: Change time','6: Change Window','7: Data Cursor','8: Keep current settings','9: End');
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
        %Go forward one volume
        elseif(options==3)
            ldx=1; ele_current=ele_save{scan_num}(ele_choice); ele_tmp=ele_save{scan_num}(ele_choice); new_scan_num=scan_num+1; ele_current=ele_save{new_scan_num}(ldx);
            while(ele_tmp~=ele_current)
                ldx=ldx+1;
                ele_current=ele_save{new_scan_num}(ldx);
                if(new_scan_num>10000)
                    fprintf('Error: no tilt exists in next volume');
                    keyboard;
                end
            end
            ele_number=ele_save{scan_num}(ele_choice);
            scan_num=new_scan_num;
            [row]=find(scan_number==scan_num);
            time_choice=file_time_save{scan_num}(1:6);
        %Go backward one volume
        elseif(options==4)
            ldx=1; ele_current=ele_save{scan_num}(ele_choice); ele_tmp=ele_save{scan_num}(ele_choice); new_scan_num=scan_num-1; ele_current=ele_save{new_scan_num}(ldx);
            while(ele_tmp~=ele_current)
                ldx=ldx+1;
                ele_current=ele_save{new_scan_num}(ldx);
                if(abs(new_scan_num)>10000)
                    fprintf('Error: no tilt exists in next volume');
                    keyboard;
                end
            end
            ele_number=ele_save{scan_num}(ele_choice);
            scan_num=new_scan_num;
            [row]=find(scan_number==scan_num);
            time_choice=file_time_save{scan_num}(1:6);
        elseif(options==5)
            fprintf(['Choose a time: \n']);
            for idx=1:max(size(file_time_save))
                fprintf([num2str(idx) ': ' num2str(roundn(ele_save(idx),-1)) ' tilt ' file_time_save{idx}(1:2) ':' file_time_save{idx}(3:4) ':' file_time_save{idx}(5:6) '\n']);
            end
            scan_num=input('Choose scan number: ');

            ele_number=ele_save{scan_num}(ele_choice);
            [row]=find(scan_number==scan_num);
            time_choice=file_time_save{scan_num}(1:6);
            
    %         if(strcmp(time_choice,'1000')
    %             if(~exist('flistd','var'))
    %                 files_tmp=dir(file_tmp);
    %                 for qdx=1:max(size(files_tmp))
    %                     tmp=files_tmp(qdx).name; tmp=char(tmp(14:19));
    %                     for ldx=1:max(size(flistd))
    %                         if(strcmp(
    %                     flistd{qdx}=char(tmp(14:19));
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
%                 if(exist('pos3','var') && big_monitor==0)
%                     figure(115)
%                     axis(lims)
%                     figure(3)
%                     axis(lims)
%                     figure(116)
%                     axis(lims)
%                     figure(4)
%                     axis(lims)
%                 else
                    figure(1)
                    axis(lims)
                    figure(2)
                    axis(lims)
                    if(big_monitor)
                        figure(3)
                        axis(lims)
                        figure(4)
                        axis(lims)
                        figure(5)
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
                    figure(5)
                    axis(lims)
                end
%             end
        end
        height_km=sqrt(range_km^2+(ae)^2+2*range_km*ae*sin(ele*pi/180))-ae+0.020;
        fprintf(['Height above ground: ' num2str(roundn(height_km,-3)) ' km' '\n']);  
    end
    