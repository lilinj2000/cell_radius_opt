function [row, col] = cellChr( freq, type )
% according the cell characteristic, group the cell
%
% the group cell is below:
%               MACRO      MICRO        None
% GSM800      (1,1)       (1,2)       (1,3)
% GSM900      (2,1)       (2,2)       (2,3)
% GSM1800     (3,1)       (3,2)       (3,3)
% GSM1900     (4,1)       (4,2)       (4,3)
% None         (5,1)       (5,2)       (5,3)

freq_list = {'GSM800', 'GSM900', 'GSM1800', 'GSM1900'};
type_list = {'MACRO', 'MICRO'};

index_freq = 5;
for ii = 1 : length(freq_list)
    if strcmpi(freq_list{ii}, freq)
        index_freq = ii;
        break;
    end
end

index_type = 3;
for ii = 1 : length(type_list)
    if strcmpi(type_list{ii}, type)
        index_type = ii;
        break;
    end
end

% assign the index to the return value
row = index_freq;
col = index_type;

end