function radius = ohPropagationAlgorithm(BSPwr, Accmin)
% The Okumura-Hata propagation model, describes the signal strength 
% decrease as a function of distance from the site. 
% The model is environment dependant as signal strength 
% from a site decreases different for urban, suburban and rural areas. 
% Furthermore antenna height is also a critical factor for propagation loss
% but is covered by assuming fixed antenna heights for each environment type. 
% Propagation loss is also dependant on frequency, 
% that is there is more loss of signal strength at 1800 MHz than at 900 MHz. 
% The above assumptions gives a signal loss of 35 dB/decade in distance.
% 
% This solution is assumed that the environment is suburban and the antenna
% height is 30 meter and the frequency is 900 MHz.

AntGain = 21;

LossdB = BSPwr + AntGain + Accmin;

radius = 1000.*power(10, (LossdB-116)./35);

end