function [nx,yt,ym,ndx]=dofault(an,fh,lp,dx,ct)

% [nx,yt,ym,ndx]=dofault(an,fh,lp,dx,ct)
% 
%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an - fault angle in degrees
% fh - fault horizontal displacement in m
% lp - length of profile in m
% dx - Min dx (the array is recalculated for power-of-two length)
% ct - crustal thickness in m
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nx - new x array
% yt - seafloor topography
% ym - mantle topography
% ndx - spacing of new x-array

%%% calculate some scalars
try     dan=deg2rad(an); % degree to radians
catch   dan=an*pi/180;
end
% fh is horizontal distance of fault, and so the vertical uplife is 
fv=fh*tan(dan); % Vertical uplift of the fault
%     
%  ^  |       /
%  |  |     t/   
%  |  |    l/
% fv  |   u/
%  |  |  a/
%  |  | f/
%  v  | /
%      ---------
%      <-- fh -->

ch=ct/tan(dan); % Horizontal distance of contact of Moho w/ fault w/respect to surface fault break
%     
%            /| ^
%          t/ | |
%         l/  | ct
%        u/   | |
%       a/    | v
%      f/<-ch->
%      /
%         
%

%%% dx is spacing, profile length is symetrical about x=0
x=-lp/2:dx:lp/2; % x array
% if lenght of array is not a power of 2, redo it
np2=power2(length(x)); % calculate next power of two length
nx=linspace(min(x),max(x),np2); % recalculate x for 2^n

%%%% make a y-array that is the same length as 'x'
% this is the topography of the crust
yt=zeros(size(nx)); % initialize y variables
% offsetting the crust topography by the crustal thickness gives mantle
% topography
ym=yt-ct;           % shift them

%%% now we find the actual intial conditions for seafloor topography
% first find all x values beyond the fault heave 
itop=find(nx>=fh/2);
% and offset them to be equal to fault vertical dispalcemnet
yt(itop)=yt(itop)+fv;
% then find the values that are on the fault surface
ifau=find(nx>-fh/2 & nx<fh/2); 
% and find their y-values according to trigonometry
yt(ifau)=[nx(ifau)-min(nx(ifau))]*tan(dan);

%%% Schematic of resultant topography 
%
%                          <--- ifau ----> <------- itop ------->
%                                         ,----------------------
%                                        /
%                                      t/   
%                                     l/
%                                    u/
%                                   a/
%                                  f/
%                                  /
%                                 /:
%                                / :
%                               /  :
%                              /   :
%                             /    :
%                            /     :
%                           /      :
%    ----------------------'       :
%                                  :
%                                 x=0                                 


%%% now do for moho topography
% we find the x's beyond the fault, which for the mantle is the fh/2, minus
% the horizontal distance we low from the crust
itop=find(nx>=fh/2-ch);
% offset this by the fault vertical displacement
ym(itop)=ym(itop)+fv;
% go find the points that lie on the fault
ifau=find(nx>-fh/2-ch & nx<fh/2-ch); 
% offset them by the trigonometry 
ym(ifau)=[nx(ifau)-min(nx(ifau))]*tan(dan)-ct;

%%% both topographies now
%                                         ,----------------------
%                                        /  crust
%                                       /------------------------   
%                                      /
%                                     /
%      sea water                     /
%                                   /    mantle
%                                  /
%                                 /:
%                                / :
%                               /  :
%                              /   :
%                             /    :
%                            /     :
%                           /      :
%    ----------------------'       :
%      crust              /        :
%    --------------------'         :
%                                  :
%                                 x=0        

%%% last is to find the new spacing after the power of 2 thing
% simply take diff of two ajacent x's
ndx=nx(2)-nx(1);

