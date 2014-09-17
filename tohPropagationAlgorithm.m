function [A, B, radius] = tohPropagationAlgorithm(cell_radius_power, toh_model)
% the trained Okumura-Hata propagation model
% 
% input:
% cell_radius_power - include the cell radius & power info
%                   - [ bspwr, accmin, real_radius]
% toh_model - in case trained, [A, B]
%
% output:
% A, B - formular power = A*log(radius) + B
% radius - the computed radius by trainded OH model

% the old OH model
% AntGain = 21;
% 
% LossdB = BSPwr + AntGain + Accmin;
% 
% radius = 1000.*power(10, (LossdB-116)./35);


if nargin<2
    toh_model = [];
end

idx_A = 1;
idx_B = 2;

idx_bspwr = 1;
idx_accmin = 2;
idx_radius = 3;

y = cell_radius_power(:, idx_bspwr) + cell_radius_power(:, idx_accmin);

if isempty(toh_model)
    

    x0 = ones(length(cell_radius_power), 1);
    x1 = log10(cell_radius_power(:, idx_radius)/1000);
    

    X = [x0, x1];
    
    [b, BINT, R, RINT, STATS] = regress(y,X);

    B = b(1);
    A = b(2);
    
    log = log4m.getLogger('toh.log');
    log.setLogLevel(log.ALL);
    log.info('tohPropagationAlgorithm()', strcat('A is ', num2str(A)));
    log.info('tohPropagationAlgorithm()', strcat('B is ', num2str(B)));
    
else
    A = toh_model(idx_A);
    B = toh_model(idx_B);
end

radius =  1000*power(10, (y-B)/A);

end