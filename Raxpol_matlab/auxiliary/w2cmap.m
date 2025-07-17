% W2CMAP WDSS-II Colormap.
%
%   W2CMAP(MAP,MIN,MAX,UNIT) generates a WDSS-II compatible colormap text
%   file using colormap MAP and the axis limits specified by MIN and MAX.
%
%   W2CMAP(MAP) generates a colormap from 0 to N-1, where N is the total
%   number of color entries in MAP.
%
%   Boon Leng Cheong
%   2/19/2009

function w2cmap(cmap,x,y,unit)
if nargin < 1
    % Demo values
    C.map = boonlib('bjetmapxinv',21);
    x = -10;
    y = 10;
    label = 'N-units';
else
    C.map = cmap;
end
if ~exist('x','var') || ~exist('y','var')
    x = 0;
    y = size(C.map,1)-1;
end
if ~exist('label','var')
    label = 'undefined';
end

C.min = x;
C.max = y;
C.unit = label;
C.num = size(C.map,1);
boonlib('cspace',C.map)

step = (C.max-C.min)/(C.num-1);

fw = fopen('wdssii_map','wt');
if (fw<0)
	fprintf('Error creating file.\n');
	return;
end

bprintf(fw,'<!--\n\n');
bprintf(fw,'  File created on %s \n',datestr(now));
bprintf(fw,'  Exported from the MATLAB script w2cmap.m \n');
bprintf(fw,'  Interested in the script? Contact boonleng@ou.edu\n\n');
bprintf(fw,'-->\n\n');

bprintf(fw,'<colorMap>\n\n');

% No data
bprintf(fw,'\t<colorBin upperBound="%.2f" name="MD">\n',C.min-10*step);
bprintf(fw,'\t<color r="0x00" g="0x00" b="0x00"></color>\n');
bprintf(fw,'\t</colorBin>\n');

% First color
c = round(C.map(1,:).'*255);
h = lower(dec2hex(c,2));
bprintf(fw,'\t<colorBin upperBound="%.2f" name="<%.2f">\n',C.min,C.min);
bprintf(fw,'\t<color r="0x%s" g="0x%s" b="0x%s"></color>\n',h(1,:),h(2,:),h(3,:));
bprintf(fw,'\t</colorBin>\n');
% Middle portion
for idx = 2:C.num-1
    u = C.min + idx*step;
    c = round(C.map(idx,:).'*255);
    h = lower(dec2hex(c,2));
    bprintf(fw,'\t<colorBin upperBound="%.2f" name="%.0f">\n',u,u-0.5*step);
    bprintf(fw,'\t<color r="0x%s" g="0x%s" b="0x%s"></color>\n',h(1,:),h(2,:),h(3,:));
    bprintf(fw,'\t</colorBin>\n');
end
% Upper portion
c = round(C.map(C.num,:).'*255);
h = lower(dec2hex(c,2));
bprintf(fw,'\t<colorBin upperBound="infinity" name=">%.2f">\n',C.max);
bprintf(fw,'\t<color r="0x%s" g="0x%s" b="0x%s"></color>\n',h(1,:),h(2,:),h(3,:));
bprintf(fw,'\t</colorBin>\n');

bprintf(fw,'\n\t<unit name="%s"></unit>\n\n',C.unit);
bprintf(fw,'</colorMap>\n');
fclose(fw);

function bprintf(varargin)
fprintf(varargin{:});
fprintf(varargin{2:end});
