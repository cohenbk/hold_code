    %Make rhohv colormap
    % x = [0 0.001 0.005 0.01 0.05 1 2 3 4 5];
    x=[0 0.1 0.2 0.3 0.4 0.5 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.92 0.94 0.96 0.98 0.99 1];
    Nx = length(x);
    clim = [min(x) max(x)];
    dx = min(diff(x));
    y = clim(1):dx:clim(2);
    for k=1:Nx-1, y(y>x(k) & y<=x(k+1)) = x(k+1); end % NEW
    cmap = colormap(boonlib('carbmap2',Nx));
    cmap2 = [...
    interp1(x(:),cmap(:,1),y(:)) ...
    interp1(x(:),cmap(:,2),y(:)) ...
    interp1(x(:),cmap(:,3),y(:)) ...
    ];