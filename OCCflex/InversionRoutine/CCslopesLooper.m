function [tslope ttopo tdistances topodistances] = CCslopesLooper(heaves,angles,crusts,infills,depths,tes)
%
% CCslopevs Te.m calculates curves of detachment slope vs distance from the axis as seen in Schouten et al  (2010) figure2')
% more te s can be added in line 6 of this m-file')
% program derived from GEOLOGYfigure2.m
% CCslopevsTe calls GEOLtestfaultnofhnote.m, flex.m, and dofault from folder FLEXURE Javier')
% Hans Schouten July 2015
% InpuTs
%   heaves: vector or scalar of horizontal extensions of the fault
%   angles: vecotr of scalar of angles for the fault geometyr
%   crusts: vector or scalar of crustal thicknesses
%   infills: vecotr or scalar of infills thicknesses
%   depths: vector or scalar of depths of fault roots
%   tes: vector of elastic thicknesses in meters
%        default = 6

% Outputs
%   slopes: matrix of slope values for Te
%   distances: matrix of coresponding y elevations
%
% loop through all Te values
% this loop generates the topography for each Te
% and it finds the slope of the topography 
% the topography is not saved, or plotted in this script.

% first initalize the counter
icount=0;

% second cd to where the calculator scripts are
cd mfiles

% now loop through parameters

for fh = heaves
    for an=angles
        for te=tes
            for ct=crusts
                for ift = infills
                    icount=icount+1;    
                    % call the calculator script
                    [fyt, ndx, nnx, yt] = slopecalc(fh,te,1,ct,an,1e5,ift);
                    % now find offset for assuming the fault begins at
                    % depth 6 km and angle = an
%                    depth=6; % km, from Schouten 
%                    dx60=depth/tan(an*pi/180);
                    % now we find the min value, and we only cut it there.
                    % this is because we had doubled the topography
                    [yy,ii]=min(fyt(1:end/2));
                    yyy=fyt(ii-4540:ii+30000);%fyt(ii:ii+30000-1);
                    xxx=ndx*(1:length(yyy));
                    % take slope of the topography
                    % not sure why we do this
                    slope=atan2(diff(yyy),diff(xxx))*180/pi;
                    % we make the slopes negative, because they are 'pointing toward' the ridge
                    % axis, That is, the footwall is sloping toward the ridge axis. 
                    tslope(:,icount)=-slope;
                    % topography of the slope data with some of the
                    % footwall
                    ttopo(:,icount)=fyt(ii-4540:ii+30000-1);
                    % full topo here
                    %ttopo2(:,icount)=fyt(1:end/2);
                    % distance vector
                    tdistances(:,icount) = nnx(1:length(tslope))'/1000;%+dx60;
                    topodistances(:,icount) = nnx(1:length(ttopo))'/1000;
                end
            end
        end
    end
end

% then cd back 

cd ..
