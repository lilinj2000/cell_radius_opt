function gsm_data = importGSMData(file)

% Header is below:
% 1~5   BSC	SITE	altitude	latitude	longitude
% 6~10  CELL	accmin	antenna_gain	antenna_tilt	antenna_type
% 11~15 bcc	bcchno	bspwr	bspwrb	c_sys_type	
% 16~20 cell_dir	cell_type	ci	env_char	height	
% 21~25 lac latitude	longitude	max_altitude	max_cell_radius	
% 26~30 mcc min_altitude	mnc	ncc	sector_angle	
% 31~32 talim override_ta

gsm_data = [];

fReadId = fopen(file);
fWriteId = fopen('GSMInput.cnai', 'w+');

header = fgetl(fReadId);
header = regexp(header, '\s', 'split');

if length(header)~=32
    % it's invalid header format
    disp('invalid header, it should be 32 column.');
    return;
end


seprator_line = fgetl(fReadId);
% count = strfind(seprator_line, '-');
if ~strcmp(seprator_line(1:3), '---')
    % invalid the seprator line
    disp('invalid seprator line, it should be begin at least 3 character ''-''');
    return;
end

while ~feof(fReadId)
    data = fgetl(fReadId);
    fprintf(fWriteId, '%s\n', data);
end

% pause(2);

fclose(fReadId);
fclose(fWriteId);

import_data = importdata('GSMInput.cnai');

gsm_data = regexp(import_data, '\s', 'split');

end