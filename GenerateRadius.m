%
% Generate Tags of normal one
% Input: measurement data in format that required by platform
% Output: radius  95% points
%         
function cell_data = GenerateRadius(measurementData, cell_data)
%function measurementData=GenerateTags(measurementData)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coloumns in MeasurementData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CGI_COL = 2;
TA_COL = 3;
LONG_COL = 17;
LAT_COL = 18;
INDEX_COL = 1;
SS_COL = [4 6 8 10 12 14 16];
NB_COL = [5 7 9 11 13 15];
TOTAL_COL_NUM = 19;
ALT_COL = 19;
TS_COL = 20;


measurementData = sortrows(measurementData, [CGI_COL]);

ClusterSchemeLayer = 1;
tagMatrix = (1:1:length(measurementData(:,1)))'*ones(1,ClusterSchemeLayer);
for i = 1:length(measurementData(:,CGI_COL))-1
    % Layer 1
    if measurementData(i, CGI_COL) == measurementData(i+1, CGI_COL)
        tagMatrix(i+1, 1) = tagMatrix(i, 1);
      
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Radius
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cgiList = unique(measurementData(:,CGI_COL));
% numOfNB=MaxGsmNeighborCellCalculated;
for numOfCgi = 1:length(cgiList(:,1))
    for Layer = 1:ClusterSchemeLayer
        CGI = cgiList(numOfCgi,:);
        CgiStartPointIndex = find(CGI==measurementData(:,CGI_COL), 1, 'first'); %Extract points of same CGI
        CgiEndPointIndex = find(CGI==measurementData(:,CGI_COL), 1, 'last');
     
        GSMCellID = CGI;
        
%         GSMCellID = GSMCellID+61106614173696;%+37937E080000
%         MCC = fix(GSMCellID/2^48)+284; %(2^56+2^52+2^51+2^50)/2^48
%         GSMCellID = GSMCellID - (MCC-284)*2^48;
%         GSMCellID = GSMCellID - fix(GSMCellID/2^44)*2^44;
%         MNC = fix(GSMCellID/2^32);
%         GSMCellID = GSMCellID-MNC * 2^32;
        % MCC=302;
        % MNC=720;
%         LAC = fix(GSMCellID/2^16);
%         CID = GSMCellID-LAC*2^16;
        
        lac = fix(GSMCellID/2^16);
        cid = GSMCellID-lac*2^16;
        
        lac_ci = strcat(num2str(lac), '-', num2str(cid)); 
        
        index_cell = 0;
        for jj=1:size(cell_data, 1)
            if strcmp(cell_data{jj, 1}, lac_ci)
                index_cell = jj;
                break;
            end
        end
        
        if index_cell==0
            continue;
        end
        
        refPoint = deg2rad(cell_data{index_cell, 2});
        
%         sizeOfCluster = CgiEndPointIndex-CgiStartPointIndex+1;
        LatLongArrayOfCluster=measurementData(CgiStartPointIndex:CgiEndPointIndex,[LAT_COL LONG_COL INDEX_COL]);
        LatLongArrayOfCluster=sortrows(LatLongArrayOfCluster,3);
        LatLongArrayOfCluster=LatLongArrayOfCluster(:,1:2);
        
        LatLongArrayOfCluster = deg2rad(LatLongArrayOfCluster);
                
        enu_lat_long = convertlatlong2enu(LatLongArrayOfCluster, refPoint);
             
        
        enu_radius = zeros(size(enu_lat_long, 1), 1);
        for jj=1:size(enu_lat_long, 1)
            enu_radius(jj) = sqrt(enu_lat_long(jj, 1).^2 + enu_lat_long(jj, 2).^2);
        end
        
        enu_radius = sort(enu_radius);
        
        cell_data{index_cell, 3} = [enu_radius(fix(size(enu_radius, 1).^0.8)), ...
            enu_radius(fix(size(enu_radius, 1).^0.95)), ...
            enu_radius(size(enu_radius, 1))];

    end
end
% if ~isempty(AllTags)
%     AllTags=sortrows(AllTags,[1 2]);
% end