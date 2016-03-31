function axis_to_use = whichaxis(axis,xmidpoint,ymidpoint,long_or_lat)
% first, find which orientation axis is
%if abs(axis(1,1)-axis(end,1)) < abs(axis(1,2)-axis(end,2))
if sum(abs(diff(axis(isfinite(axis(:,1)),1)))) < sum(abs(diff(axis(isfinite(axis(:,2)),2))))
    MA = 1;
    MM = 2;
else
    MA = 2;
    MM = 1;
end

[nul mid2ix] = min(abs(long_or_lat - ymidpoint));

% do axes stuff
axsc = 1;
ax1 = NaN*zeros(size(axis,1),size(axis,2)/2)';
ax2 = ax1;
for jk = 1:2:size(axis,2)
    ax1(axsc,find(isfinite(axis(:,jk-1+MA)))) = axis(find(isfinite(axis(:,jk-1+MA))),jk-1+MA);
    ax2(axsc,find(isfinite(axis(:,jk-1+MM)))) = axis(find(isfinite(axis(:,jk-1+MM))),jk-1+MM);
    ax1midvals(axsc) = ax1(axsc,mid2ix);
    axsc = axsc+1;
end

[nul ax1idx ] = min(abs(xmidpoint - ax1midvals));
axis_to_use(:,1) = ax1(ax1idx,:)';
axis_to_use(:,2) = ax2(ax1idx,:)';