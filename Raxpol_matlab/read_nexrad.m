%Read NEXRAD Data

cd(dir_name)

%Get global attributes
ncid=files(row).name;
data.Station=ncreadatt(ncid,'/','instrument_name');
data.time_cov_start=ncreadatt(ncid,'/','time_coverage_start');
data.time_cov_end=ncread(ncid,'time_coverage_end');
data.sweep_num=ncread(ncid,'sweep_number');
data.sweep_mode=ncread(ncid,'sweep_mode');
data.fixed_angle=ncread(ncid,'fixed_angle');
data.time=ncread(ncid,'time');
data.range=ncread(ncid,'range');
data.az=ncread(ncid,'azimuth');
data.el=ncread(ncid,'elevation');
data.pw=ncread(ncid,'pulse_width');
data.prt=ncread(ncid,'prt');
data.va=nanmean(ncread(ncid,'nyquist_velocity'));
data.lat=nanmean(ncread(ncid,'latitude'));
data.lon=nanmean(ncread(ncid,'longitude'));
%data.heading=nanmean(ncread(ncid,'heading'));
data.Z_H=ncread(ncid,'REF');
data.Z_DR=ncread(ncid,'ZDR');
data.vr_h=ncread(ncid,'VEL');
data.Phi_DP=ncread(ncid,'PHI');
data.rho_hv=ncread(ncid,'RHO');
data.SW=ncread(ncid,'SW');
%data.Z_corr=ncread(ncid,'DBZC');
%data.ZDR_corr=ncread(ncid,'ZDRC');

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