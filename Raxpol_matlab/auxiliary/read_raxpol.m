% %Read RaXPol Data
% 
% cd(dir_name)
% %cd(radar_path)
% 
% %Get global attributes
% ncid=files(row).name;
% %ncid=radar_file;
% data.time_cov_start=ncread(ncid,'time_coverage_start');
% data.time_cov_end=ncread(ncid,'time_coverage_end');
% data.sweep_num=ncread(ncid,'sweep_number');
% data.sweep_mode=ncread(ncid,'sweep_mode');
% data.fixed_angle=ncread(ncid,'fixed_angle');
% %data.altitude_corr=ncread(ncid,'altitude_correction');
% data.altitude=ncread(ncid,'altitude');
% data.altitude_agl=ncread(ncid,'altitude_agl');
% data.time=ncread(ncid,'time');
% data.range=ncread(ncid,'range');
% data.az=ncread(ncid,'azimuth');
% data.el=ncread(ncid,'elevation');
% data.pw=ncread(ncid,'pulse_width');
% data.prt=ncread(ncid,'prt');
% data.va=nanmean(ncread(ncid,'nyquist_velocity'));
% data.lat=nanmean(ncread(ncid,'latitude'));
% data.lon=nanmean(ncread(ncid,'longitude'));
% % data.heading=nanmean(ncread(ncid,'heading'));
% data.Z_H=ncread(ncid,'DBZ');
% % data.Z_DR=ncread(ncid,'ZDR');
% 
% % data.Z_H=ncread(ncid,'REF');
% data.Z_DR=ncread(ncid,'ZDR');
% data.vr_h=ncread(ncid,'VEL');
% %%% Dodge City uses VUN created by B. Cohen
% % data.vr_h=ncread(ncid,'VUN');
% %%% El Reno uses UV created by J. Houser
% % data.vr_h=ncread(ncid,'VU');
% %data.vr_hc=ncread(ncid,'VC');
% %data.SW_h=ncread(ncid,'WIDTH');
% % data.Phi_DP=ncread(ncid,'PHIDP');
% % data.Phi_DP=ncread(ncid,'PHI'); %WSR
% % data.KDP=ncread(ncid,'KDP');
% data.rho_hv=ncread(ncid,'RHOHV'); %'RHOHV'
% %data.Z_corr=ncread(ncid,'ZC');
% %data.ZDR_corr=ncread(ncid,'DRC');
% 
% data.detr=data.range(2)-data.range(1);
% data.gw=data.detr;
% 
% 
% % data.Z_H(data.Z_H==data.missing_data)=NaN;
% % data.Z_DR(data.Z_DR==data.missing_data)=NaN;
% % data.vr_h(data.vr_h==data.missing_data)=NaN;
% % data.SW(data.SW==data.missing_data)=NaN;
% % data.rho_hv(data.rho_hv==data.missing_data)=NaN;
% % data.phi_dp(data.phi_dp==data.missing_data)=NaN;
% 
% vr_lims=[-nanmax(abs(data.vr_h(:))) nanmax(abs(data.vr_h(:)))];
% 
% cd(base_dir)

%Read Horus Data

cd(dir_name)
%cd(radar_path)

%Get global attributes
ncid=files(row).name;
%ncid=radar_file;
data.time_cov_start=ncread(ncid,'time_coverage_start');
data.time_cov_end=ncread(ncid,'time_coverage_end');
data.sweep_num=ncread(ncid,'sweep_number');
data.sweep_mode=ncread(ncid,'sweep_mode');
data.fixed_angle=ncread(ncid,'fixed_angle');
%data.altitude_corr=ncread(ncid,'altitude_correction');
data.altitude=ncread(ncid,'altitude');
data.altitude_agl=ncread(ncid,'altitude_agl');
data.time=ncread(ncid,'time');
data.range=ncread(ncid,'range');
data.az=ncread(ncid,'azimuth');
data.el=ncread(ncid,'elevation');
data.pw=ncread(ncid,'pulse_width');
data.prt=ncread(ncid,'prt');
data.va=nanmean(ncread(ncid,'nyquist_velocity'));
data.lat=nanmean(ncread(ncid,'latitude'));
data.lon=nanmean(ncread(ncid,'longitude'));
% data.heading=nanmean(ncread(ncid,'heading'));
data.Z_H=ncread(ncid,'reflecti');
% data.Z_DR=ncread(ncid,'ZDR');

% data.Z_H=ncread(ncid,'REF');
data.Z_DR=ncread(ncid,'differen');
data.vr_h=ncread(ncid,'VUN2');
% data.vr_h=ncread(ncid,'V_u2');
%%% Dodge City uses VUN created by B. Cohen
% data.vr_h=ncread(ncid,'VUN');
%%% El Reno uses UV created by J. Houser
% data.vr_h=ncread(ncid,'VU');
%data.vr_hc=ncread(ncid,'VC');
%data.SW_h=ncread(ncid,'WIDTH');
% data.Phi_DP=ncread(ncid,'PHIDP');
% data.Phi_DP=ncread(ncid,'PHI'); %WSR
% data.KDP=ncread(ncid,'KDP');
data.rho_hv=ncread(ncid,'cross_co'); %'RHOHV'
%data.Z_corr=ncread(ncid,'ZC');
%data.ZDR_corr=ncread(ncid,'DRC');

data.detr=data.range(2)-data.range(1);
data.gw=data.detr;


% data.Z_H(data.Z_H==data.missing_data)=NaN;
% data.Z_DR(data.Z_DR==data.missing_data)=NaN;
% data.vr_h(data.vr_h==data.missing_data)=NaN;
% data.SW(data.SW==data.missing_data)=NaN;
% data.rho_hv(data.rho_hv==data.missing_data)=NaN;
% data.phi_dp(data.phi_dp==data.missing_data)=NaN;

vr_lims=[-nanmax(abs(data.vr_h(:))) nanmax(abs(data.vr_h(:)))];

cd(base_dir)