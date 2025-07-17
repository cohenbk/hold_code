dnum=datenum(data.time_cov_start','yyyy-mm-ddTHH:MM:SSZ');
if(ele_save(scan_num)>=10)
    fname_save=[data.Station '_' datestr(dnum,'yyyymmdd') '_' datestr(dnum,'HHMMSS') '_' num2str(round(ele_save(scan_num)*10)/10)];
else
    fname_save=[data.Station '_' datestr(dnum,'yyyymmdd') '_' datestr(dnum,'HHMMSS') '_' strcat('0',num2str(round(ele_save(scan_num)*10)/10))]; 
end
eval('cd stats')

%Added code to get user radius
while(rthres==-1)
    rthres=input('Enter desired radius: ');
    fprintf(['The radius you have selected is: ' num2str(rthres) '\n'])
    accept_radius=input('Enter 1 to accept, 0 to redo: ');
    if(~accept_radius)
        rthres=-1;
    end
end
if(~exist(fname_save,'file'))
    pause;
    figure(2)
    [xs,ys]=ginput(1);
    rr=sqrt((xx-xs).^2+(yy-ys).^2);
    yn=rr<rthres;
    rho_save=data.rho_hv(yn);
    Zh_tmp=data.Z_H(1:size(data.vr_h,1),1:size(data.vr_h,2));
    Z_save=Zh_tmp(yn);
    vr_save=data.vr_h(yn);
    clear Zh_tmp
    Zdr_save=data.Z_DR(yn);
    data_el=nanmean(data.el);
    save(fname_save,'rho_save','Z_save','Zdr_save','vr_save','rthres','xs','ys','data_el');
else
    file_ov=input('Overwrite existing file? ');
    if(file_ov)
        figure(2)
        [xs,ys]=ginput(1);
        if(~exist('rthres','var'))
            rthres=input('Choose radius');
        end
        rr=sqrt((xx-xs).^2+(yy-ys).^2);
        yn=rr<rthres;
        rho_save=data.rho_hv(yn);
        Zh_tmp=data.Z_H(1:size(data.vr_h,1),1:size(data.vr_h,2));
        Z_save=Zh_tmp(yn);
        vr_save=data.vr_h(yn);
        clear Zh_tmp
        Zdr_save=data.Z_DR(yn);
        data_el=nanmean(data.el);
        %Save file with desired data within radius
        save(fname_save,'vr_save','rho_save','Z_save','Zdr_save','rthres','xs','ys','data_el');
    end
end
eval('cd ..');
options=8;