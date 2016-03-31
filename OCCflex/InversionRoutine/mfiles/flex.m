function w=flex(ht,hc,dx,te)
% gives response to load in meters (- down)
% w=flex(ht,hc,dx,te)
%%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ht=water-crust interface topo, centered at 0, m (array)
% hc=crust-mantle interface topo, centered at 0, m (array)
% dx=spacing in m
% te= Elastic thickness, in m
%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w - the response of the plate
% assumes:
% E=1e11;  Young's modulus
% v=.25 Poisson's ratio
% rho water     = 1.03 
% rho crust     = 2.7
% rho mantle    = 3.3
% weissel & karner 1989, jgr 94 13919-13950

%E=Emult*1e10; % mlarson change, aug 20 to test different E
E=1e11;                 % Youngs modulus
v=.25;                  % Poisson's ratio
g=9.81;                 % Gravity
rho_w=1030;             % water
rho_c=2700;             % crust
rho_m=3300;             % mantle density
rho_s=3000;             % serpentinite density

% eqn 3.115 from Turcotte & Schubert 3rd Ed.
% isostatic result for ht and hc
% this is units of elevation
st=ht*(rho_c-rho_w)/(rho_m-rho_w);
sc=hc*(rho_m-rho_c)/(rho_m-rho_w);
% these results look like this 
%
%
%                      ,------------,
%                     ,              ,
%                    /                \
%                  ,                    ,
%                 /                      \
%               ,      ,------------,      ,
%              /    ,'                ',    \ 
% ------------' , '                      ',  '------------ st
% -----------'                                '------------ sc

% now add and invert
s=-(sc+st);
% the result is a steeper slope than before
%
% ------------,                          ,------------ s
%              \                        /
%               \                      /
%                \                    /
%                 \                  /
%                  \                /
%                   \              /
%                    '------------'

% find some wavelenghts, and their spectrum
[S,k]=jfft(s,dx);
% x = 100000*[1:length(k)]/length(k);
% E = E*(1+x/10);

% find D, a scalar, whic is the rigidity of the plate. 
% from Turcotte & Schubert 3rd Ed. eqn. 3.72, also, Buck 1988
D=E.*te^3/(12*(1-v^2));

% find compensation, from Turcotte Schubert 3rd Ed. eqn 3.117
cte=(rho_m-rho_w)*(D/g.*k.^4+((rho_m-rho_w))).^(-1);
% this is a curve that starts at x=0, y = 1, 
% and then goes to y=0 as x -> length profile
% curvature is changed with all the factors, but Te is the one we deal with
% -.
%   `
%   :
%   :
%   :
%   :
%    `---------------------------------------------------------------- 
% 0                                                         

W=cte.*S';
% now time to do the inverse fft, and find the real part
% now we have our w
w=jifft(W);
