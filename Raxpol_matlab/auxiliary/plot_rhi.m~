continue_rhi_flag=true; rhi_fig=false;
while(continue_rhi_flag)
    
    if(~exist('input_azi','var'))
        input_azi=input('Choose azimuth angle: ');
        input_storm_motion=input('Type storm motion (dr/dt): ');
        start_range=input('Starting range (km): ');
        end_range=input('Ending range (km): ');
    end

    cd(file_dir)
        files_tmp=dir(['RAXPOL_' date_save{scan_num} '-' time_choice '*']);
    cd(root_dir)
    max_eles=max(size(files_tmp));
    for rdx=1:max_eles
        if(data_flags(2)==1)
            if(rdx<11)
                filelist=['RAXPOL_' date_save{scan_num} '-' time_choice '_0' num2str(rdx-1) '.mat'];
            else
                filelist=['RAXPOL_' date_save{scan_num} '-' time_choice '_' num2str(rdx-1) '.mat'];
            end
            fprintf(['Loading : ' filelist '\n'])
        end
        cd(file_dir)
            load(filelist);
        cd(root_dir)
        if(minus_sign_vr)
            data.vr_h=-data.vr_h;
        end
        [c1,azi_choice]=min(abs(data.az-input_azi));
        [c1,r1]=min(abs(start_range-r));
        [c2,r2]=min(abs(end_range-r));
        if(isfield('data','Z_corr'))
            Z_DR(1:r2-r1+1,rdx)=data.Zdr_corr(r1:r2,azi_choice);
            Z_H(1:r2-r1+1,rdx)=data.Z_corr(r1:r2,azi_choice);
        else   
            Z_DR(1:r2-r1+1,rdx)=data.Z_DR(r1:r2,azi_choice);
            Z_H(1:r2-r1+1,rdx)=data.Z_H(r1:r2,azi_choice);
        end
        vr(1:r2-r1+1,rdx)=data.vr_h(r1:r2,azi_choice);
        rho_hv(1:r2-r1+1,rdx)=data.rho_hv(r1:r2,azi_choice);
        time_diff=datenum([0 0 0 str2num(file_time_save{scan_num+1}(1:2)) str2num(file_time_save{scan_num+1}(3:4)) str2num(file_time_save{scan_num+1}(5:6))])-datenum([0 0 0 str2num(file_time_save{scan_num}(1:2)) str2num(file_time_save{scan_num}(3:4)) str2num(file_time_save{scan_num}(5:6))]); 
        time_diff=time_diff/max_eles;
        dist_shift(rdx)=input_storm_motion*time_diff*(rdx-1);
        ele_rhi(rdx)=nanmean(data.el);
    end    

    rr=repmat(r',1,max_eles);
    ee=repmat(ele_rhi,max(size(r)),1);
    dist_shift=repmat(dist_shift,max(size(r)),1);
    rr=rr+dist_shift;
    zz=sqrt(rr.^2+(ae)^2+2*rr.*ae.*sin(ee.*pi/180))-ae;
    if(~rhi_fig)
        figure(10)
        feval('boonlib','bsizewin',10,winsize)
%         hc1=subplot(3,1,1)
        subplot(3,1,1)
        data_tmp=double(Z_H);
        hr1=pcolor(rr(r1:r2,:),zz(r1:r2,:),data_tmp);
        colormap(carbone42)
        caxis(clims3)
        title(['^{\circ} Reflectivity (dBZ) at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','tex','HorizontalAlignment','Center')
        xlabel('Zonal Distance (km)','FontSize',axis_label_size)
        ylabel('Meridonal Distance (km)','FontSize',axis_label_size)
        colorbar
        % set(gca,'DataAspectRatio',[1 1 1])
        shading flat
        subplot(3,1,2)
        data_tmp=double(Z_DR);
        hr2=pcolor(rr(r1:r2,:),zz(r1:r2,:),data_tmp);
        colormap(carbone42)
        caxis([zdr_lims])
        title(['^{\circ} Z_{DR} at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','tex','HorizontalAlignment','Center')
        xlabel('Zonal Distance (km)','FontSize',axis_label_size)
        ylabel('Meridonal Distance (km)','FontSize',axis_label_size)
        colorbar
        % set(gca,'DataAspectRatio',[1 1 1])
        shading flat
        subplot(3,1,3)
        data_tmp=double(rho_hv);
        hr3=pcolor(rr(r1:r2,:),zz(r1:r2,:),data_tmp);
        colormap(carbone42)
        caxis([0.9 1])
        title(['\rho_{hv} at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','tex','HorizontalAlignment','Center')
        xlabel('Zonal Distance (km)','FontSize',axis_label_size)
        ylabel('Meridonal Distance (km)','FontSize',axis_label_size)
        colorbar
        % set(gca,'DataAspectRatio',[1 1 1])
        shading flat
        rhi_fig=true;
    else
        figure(10)
        subplot(3,1,1)
        set(hr1,'XData',rr(r1:r2,:))
        set(hr1,'YData',zz(r1:r2,:))
        set(hr1,'CData',double(Z_H))
        subplot(3,1,2)
        set(hr2,'XData',rr(r1:r2,:))
        set(hr2,'YData',zz(r1:r2,:))
        set(hr2,'CData',double(Z_DR))
        subplot(3,1,3)
        set(hr3,'XData',rr(r1:r2,:))
        set(hr3,'YData',zz(r1:r2,:))
        set(hr3,'CData',double(rho_hv))
    end

    clear rr ee zz c1 c2 i1 i2 Z_H Z_DR rho_hv vr ele_rhi dist_shift;
%     rhi_command=input('(1) +1 volume, (2) -1 volume, (0) stop');
    figure(103)
    set(103,'Position',pos9)
    h11 = uicontrol('Style', 'pushbutton', 'String', 'Go forward 1 vol','Position', [2 2 100 50], 'Callback', 'rhi_command=1; uiresume');
    h21 = uicontrol('Style', 'pushbutton', 'String', 'Go back 1 vol','Position', [102 2 100 50], 'Callback', 'rhi_command=2; uiresume');
    h31 = uicontrol('Style', 'pushbutton', 'String', 'Stop RHIs','Position', [202 2 100 50], 'Callback', 'rhi_command=0; uiresume');
    uiwait(103) 
    if(rhi_command==1)
        scan_num=scan_num+1;
        [row]=find(scan_number==scan_num);
        time_choice=file_time_save{scan_num}(1:6);
    elseif(rhi_command==2)
        scan_num=scan_num-1;
        [row]=find(scan_number==scan_num);
        time_choice=file_time_save{scan_num}(1:6);
    else
        continue_rhi_flag=false;
    end
end
clear input_azi input_storm_motion start_range end_range rhi_fig hr1 hr2 hr3; close(103);