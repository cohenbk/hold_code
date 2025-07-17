function idx = raz2idx(dr,nr,az_set,tol,x,y)

[xx yy] = meshgrid(x, y);

az = atan2(xx,yy)/pi*180;
az(az<0) = az(az<0)+360;

range = sqrt(xx.^2+yy.^2);

az_set(az_set<0) = az_set(az_set<0)+360;

ia = zeros(size(az));
ma = ia;
for k = 1:numel(az)
	tmp = az_set-az(k);
	[ma(k) ia(k)] = min(abs(tmp));
end
ir = floor(range/dr)+1;
idx = (ia-1)*nr+ir;
idx(ma>tol) = 1;
