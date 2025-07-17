function [x,y]=latlon2cart(latcoords,loncoords,latradar,lonradar)

r_earth=6378.1;
phi_s=latcoords*pi/180; lambda_s=loncoords*pi/180; phi_f=latradar*pi/180; lambda_f=lonradar*pi/180; 
% x=repmat(cos(phi_f)*cos(lambda_f),1,max(size(latcoords)))-cos(phi_s).*cos(lambda_s);
% y=repmat(cos(phi_f)*sin(lambda_f),1,max(size(loncoords)))-cos(phi_s).*sin(lambda_s);
dlat=phi_f-phi_s;
dlon=lambda_f-lambda_s;
a=sin(dlat/2).*sin(dlat/2)+(cos(phi_s)*cos(phi_f)).*(sin(dlon/2).*sin(dlon/2));
c=2*atan2(sqrt(a),sqrt(1-a));
d=r_earth*c;
y1=sin(dlon)*cos(phi_f);
x1=cos(phi_s)*sin(phi_f)-sin(phi_s)*cos(phi_f).*cos(dlon);
theta=atan2(y1,x1);

%Convert to radar coordinates
x=d.*sin(theta+pi);
y=d.*cos(theta+pi);

end

%OU-PRIME: 35.1801,-97.4337