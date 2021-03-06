function deg = dms2deg( dms )
% dms2deg just convert the dms(degree, minus, second) to degree
% the format maybe EDDMMSS, WDDMMSS, SDDMMSS, NDDMMSS;
% DDMMSSE, DDMMSS

if iscell(dms)
    dms = char(dms);
end

if dms(1)=='N' ...
    || dms(1)=='S' ...
    || dms(1)=='E' ...
    || dms(1)=='W'

    sig = dms(1);
    dms(1) = [];
    
    dms(length(dms)+1) = sig;
    deg = str2angle(dms);
else
    deg = str2double(dms);
end




end