clear;
clc;

GSMImport = importdata('GSM.cnai');

% GSMImport = importdata('all.cnai');

% GSM = cell(length(GSMImport), 1);
for ii = 1: length(GSMImport)
    GSM{ii} = regexp(GSMImport{ii}, '\s', 'split');
end

GSM = GSM';


% Header = GSM{1};
% Header = regexp(Header, '\s', 'split');
% 
% Header(2)
% 
% Header(3)