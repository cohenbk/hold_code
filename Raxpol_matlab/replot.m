cd(dir_name)
    read_raxpol;
cd(base_dir)

if(sum(ismember(plot_choice,9))>0)
    data.SW=data.SW_h;
    sw_lims=[0 15];
end
if(KOUN_adjust_flag)
   zdr_lims=[-5 8]; 
end
clims3=[-15 75];

%     if(radar_choice<3)
%Find imaginary values and set to NaN
[row2]=find(abs(imag(data.rho_hv))>0);
data.rho_hv(row2)=NaN;
if(att_corr_flag)
[row2]=find(abs(imag(data.ZDR_corr))>0);
data.ZDR_corr(row2)=NaN; clear row2;
data.ZDR_corr=data.ZDR_corr+zdr_calib;
data.ZDR_corr(data.ZDR_corr<-900)=NaN;
data.Z_corr(data.Z_corr<-900)=NaN;
end
data.Z_H(data.Z_H<-900)=NaN;
data.Z_DR(data.Z_DR<-900)=NaN;
if(isfield(data,'phi_dp'))
    data.phi_dp(data.phi_dp<-900)=NaN;
end

if(isfield(data,'phi_dp_filtered'))
    data.phi_dp_filtered(data.phi_dp_filtered<-900)=NaN;
end
data.rho_hv(data.rho_hv<-900)=NaN;
%     end

% %Some things to add if it is RaXPol data
% if(isfield(data,'roll'))
%     data.detr=data.range(2)-data.range(1);
%     if(data.heading<85)
%         data.heading=data.heading+8.7;
%     else
%         data.heading=data.heading-8.7;
%     end
% %     if(strcmp(data.time_start(1:10)','2011-06-19'))
% %         data.heading=90.8236; %Override heading due to GPS problem
% %     end
% end

if(isfield(data,'va'))
    vr_lims=[-data.va data.va];
end

%If a radar name hasn't been provided in the file, it will ask for it
%here

if(min(size(strfind(dir_name,'SPol')))>0)
    data.Station='SPol'; 
end
if(~isfield(data,'Station'))
    if(isfield(data,'roll') || isfield(data,'heading'))
        data.Station='RaXPol';
    else
        data.Station=input('Input radar name: ');
    end
end
radar=data.Station;
%     keyboard;

az=data.az; 
if(isfield(data,'el'))
    ele=median(data.el);
else
    ele=median(data.ele);
end
%     
%     if(size(data.rho_hv,2)==max(size(az)))
%         num_gates=size(data.rho_hv,1);
%     else
%         num_gates=size(data.rho_hv,2);
%     end

if(size(data.vr_h,2)==max(size(az)))
    num_gates=size(data.vr_h,1);
else
    num_gates=size(data.vr_h,2);
end

if(strcmp(data.Station,'SPol'))
    r=data.r/1000; 
else
    if(ouprime_flag || strcmp(data.Station,'RaXPol'))
        r=data.detr+data.detr*([1:num_gates]-1); r=r/1000;
    else
        r=data.detr*8+data.detr*([1:num_gates]-1); r=r/1000;
    end
end

% if(isfield(data,'heading'))
% %     if(data.heading>90) 
% %         data.heading=0;
% %         %data.heading=82.1017;
% %     end
%     az=az-data.heading;
%     az(az>360)=az(az>360)-360;
%     az(az<0)=az(az<0)+360;
% end
az_rad=az*pi/180;
az_rad=repmat(az_rad,1,max(size(r)));

r_km=repmat(r,size(data.az,1),1);
if(data_cursor_mode)
    xx=r_km.*sin(az_rad)*cos(ele*pi/180); yy=r_km.*cos(az_rad)*cos(ele*pi/180);
    if(isfield(data,'az2'))
        az_rad=data.az2*pi/180;
        az_rad=repmat(az_rad,1,max(size(r)));
        r_km=repmat(r,size(data.az2,1),1);
        xx2=r_km.*sin(az_rad); yy2=r_km.*cos(az_rad);
    elseif(isfield(data,'az_vr'))
        az_rad=data.az_vr*pi/180;
        az_rad=repmat(az_rad,1,max(size(r)));
        r_km=repmat(r,size(data.az_vr,1),1);
        xx2=r_km.*sin(az_rad)*cos(ele*pi/180); yy2=r_km.*cos(az_rad)*cos(ele*pi/180);
        xx2=xx2'; yy2=yy2';
    end
else
    xx=r_km.*sin(az_rad)*cos(ele*pi/180); yy=r_km.*cos(az_rad)*cos(ele*pi/180);
%         if(~ouprime_flag)
%              xx=xx'; yy=yy';
%         end
    if(isfield(data,'az2'))
        az_rad=data.az2*pi/180;
        az_rad=repmat(az_rad,1,max(size(r)));
        r_km=repmat(r,size(data.az2,1),1);
        xx2=r_km.*sin(az_rad)*cos(ele*pi/180); yy2=r_km.*cos(az_rad)*cos(ele*pi/180);
    elseif(isfield(data,'az_vr'))
        az_rad=data.az_vr*pi/180;
        az_rad=repmat(az_rad,1,max(size(r)));
        r_km=repmat(r,size(data.az2,1),1);
        xx2=r_km.*sin(az_rad)*cos(ele*pi/180); yy2=r_km.*cos(az_rad)*cos(ele*pi/180);
    end
end
    %Transpose data if needed
    if(size(xx,1)~=size(data.vr_h,1))
        xx=xx'; yy=yy';
    end
%         if(size(xx2,1)~=size(data.vr_h,1))
%             xx2=xx2'; yy2=yy2';
%         end


%     if(KOUN_adjust_flag && radar_choice==1)
%        xx=xx+2.6; yy=yy-6.2;
%        if(exist('xx2','var'))
%             xx2=xx2+2.6; yy2=yy2-6.2;
%        end
%     end
%     [data.Z_H,data.Z_DR]=att_corr(data.Z_H,data.Z_DR,data.phi_dp,data.rho_hv);
num_loops=max(size(plot_choice));

for pdx=1:num_loops
    figure(pdx)
    if(plot_choice(pdx)==7 || (plot_choice(pdx)>=9 && exist('xx2','var')) || (plot_choice(pdx)==2 && exist('xx2','var')))
        set(p1(pdx),'XData',cat(2,xx2,xx2(:,1)));
        set(p1(pdx),'YData',cat(2,yy2,yy2(:,1)));
    else
        set(p1(pdx),'XData',cat(2,xx,xx(:,1)));
        set(p1(pdx),'YData',cat(2,yy,yy(:,1)));
    end
    xtmp=cat(2,xx,xx(:,1)); ytmp=cat(2,yy,yy(:,1));
    if(exist('xx2','var'))
        xtmp2=cat(2,xx2,xx2(:,1)); ytmp2=cat(2,yy2,yy2(:,1));
    end
    set(p1(pdx),'ZData',double(zeros(size(xtmp,1),size(xtmp,2))))
    if(plot_choice(pdx)==1)
        if(att_corr_flag)
            data_tmp=double(data.Z_corr);
            set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)))
        else
            data_tmp=double(data.Z_H(1:size(data.vr_h,1),1:size(data.vr_h,2)));
            set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)))
        end
        t1(pdx)=title('Hi');
        set(t1(pdx),'String',[num2str(round(ele*10)/10) '$^{\circ}$ Z (dBZ) at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','Latex','HorizontalAlignment','Center')                
        if(new_color_flag(pdx))
            colormap(cmap3)
            caxis(clims3)
            colorbar
            new_color_flag(pdx)=false;
        end
    elseif(plot_choice(pdx)==2)
        if(isfield(data,'vr_dealiased'))
            data_tmp=double(data.vr_dealiased);
            set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1))); 
        elseif(isfield(data,'fold_int_map'))
            data_tmp=double(data.vr_h)+double(data.fold_int_map*data.va);
            set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1))); 
        else
            data_tmp=double(data.vr_h);
            set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)));
        end

        t1(pdx)=title('Hi');
        set(t1(pdx),'String',[num2str(round(ele*10)/10) '$^{\circ}$ $v_{r}$ at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','LaTex','HorizontalAlignment','Center')                
        if(new_color_flag(pdx))
            colormap(cmap5)
            if(isfield(data,'vr_dealiased') || isfield(data,'fold_int_map'))
                caxis(vr_lims2)
            else
                caxis(vr_lims2)
            end
            colorbar
            new_color_flag(pdx)=false;
        end
    elseif(plot_choice(pdx)==3)
        data_tmp=double(data.rho_hv);
        set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)))
        t1(pdx)=title('Hi');
        set(t1(pdx),'String',[num2str(round(ele*10)/10) '$^{\circ}$ $\rho_{  HV}$ at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','Latex','HorizontalAlignment','Center')                            
        if(new_color_flag(pdx))
            colormap(cmap2)
            caxis([0 1])
            colorbar
            new_color_flag(pdx)=false;
        end
    elseif(plot_choice(pdx)==4)
        if(att_corr_flag)
            data_tmp=double(data.ZDR_corr);
            set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)));
        else
            data_tmp=double(data.Z_DR);
            set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)));
        end
        t1(pdx)=title('Hi');
        set(t1(pdx),'String',[num2str(round(ele*10)/10) '$^{\circ}$ $Z_{DR}$ at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','Latex','HorizontalAlignment','Center')                
        if(new_color_flag(pdx))
            colormap(cmap3)
            caxis(zdr_lims)
            colorbar
            new_color_flag(pdx)=false;
        end
    elseif(plot_choice(pdx)==5)
        if(att_corr_flag)
            yn1=data.Z_corr>Z_thres&data.rho_hv<rhohv_thres; yn2=data.Z_corr>Z_thres&data.rho_hv<0.5; yn3=data.Z_corr>Z_thres&data.rho_hv<0.3; yn4=data.Z_corr>Z_thres&data.rho_hv<0.1;
        else
            yn1=data.Z_H(1:size(data.vr_h,1),1:size(data.vr_h,2))>Z_thres&data.rho_hv<rhohv_thres; yn2=data.Z_H(1:size(data.vr_h,1),1:size(data.vr_h,2))>Z_thres&data.rho_hv<0.5; yn3=data.Z_H(1:size(data.vr_h,1),1:size(data.vr_h,2))>Z_thres&data.rho_hv<0.3; yn4=data.Z_H(1:size(data.vr_h,1),1:size(data.vr_h,2))>Z_thres&data.rho_hv<0.1;
        end
        data_tmp=double(yn1)+double(yn2)+double(yn3)+double(yn4); 
        set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)))
        title([num2str(round(ele*10)/10) '^{\circ} Thresholded image at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','tex','HorizontalAlignment','Center')
        if(new_color_flag(pdx))
            caxis([0 4])
            colormap(jet(5))
            hcb=colorbar('YTickLabel',...
            {'>0.82';'0.5-0.82';'0.3-0.5';'0.1-0.3';'<0.1'});
            set(hcb,'YTick',0.4:0.8:3.6)
            new_color_flag(pdx)=false;
        end
    elseif(plot_choice(pdx)==6)
        if(att_corr_flag)
            data_tmp=double(data.phi_dp_filered);
            set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)))
        else
            data_tmp=double(data.phi_dp);
            set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)))
        end
        t1(pdx)=title('Hi');
        set(t1(pdx),'String',[num2str(round(ele*10)/10) '^{\circ} \Phi_{DP} at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','tex','HorizontalAlignment','Center')                
        colormap(cmap3)
        caxis([-180 180])
    elseif(plot_choice(pdx)==7)
        data_tmp=double(data.KDP);
        set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)))
        title([num2str(round(ele*10)/10) '^{\circ} KDP at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','tex','HorizontalAlignment','Center')
        colormap(cmap3)
        caxis([-2 10])
    elseif(plot_choice(pdx)==8)
        data_tmp=double(data.A_H);
        set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)))
        title([num2str(round(ele*10)/10) '^{\circ} A_{H} at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','tex','HorizontalAlignment','Center')
        colormap(cmap3)
        caxis([-2 15])
    elseif(plot_choice(pdx)==9)
        data_tmp=double(data.SW);
        set(p1(pdx),'CData',cat(2,data_tmp,data_tmp(:,1)))
        title([num2str(round(ele*10)/10) '^{\circ} SW at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','tex','HorizontalAlignment','Center')
        colormap(cmap3)
    %                 caxis([-2 15])
        caxis(sw_lims)
    elseif(plot_choice(pdx)==10)
        set(p1(pdx),'CData',double(sd_phidp(1:size(data.vr_h,1),1:size(data.vr_h,2))))
        title([num2str(round(ele*10)/10) '^{\circ} \sigma_{\Phi_{DP}} at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','tex','HorizontalAlignment','Center')
        colormap(cmap3)
    elseif(plot_choice(pdx)==11)
        set(p1(pdx),'CData',double(sd_ZDR(1:size(data.vr_h,1),1:size(data.vr_h,2))))
        title([num2str(round(ele*10)/10) '^{\circ} \sigma_{\Z_{DR}} at ' time_choice ' UTC'],'FontSize',title_size,'Interpreter','tex','HorizontalAlignment','Center')
        colormap(cmap3)
    end
end
