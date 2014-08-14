function [inner_radius, outer_radius] = taAlgorithm( ta )
% the TA value holds information about the delay in the air  interface
% between MS and the base station.
%  1 TA ~ 550 m
% 
% TA Look-Up Table
% TA  Inner radius    Outer radius      diff    TA
% 0     0               450
% 1     0               825             825     1.5
% 2     275             1375            1100    2
% 3     825             1925            1100    2
% 4     1375            2475            1100    2
% ...
% 63    33825           34925           1100    2

taLimStep = 120000/219; % 547.9425

if ta>=2
    inner_radius = taLimStep * (ta-1.5);
    outer_radius = taLimStep * (ta+0.5);
elseif ta==1
    inner_radius = 0;
    outer_radius = 825;
elseif ta==0
    inner_radius = 0;
    outer_radius = 450;   
end

end