function drawCellMsrGeo(cell_info, cell_angle)

CGI_COL = 2;
LONG_COL = 17;
LAT_COL = 18;
INDEX_COL = 1;

% cgiList = unique(measurement(:,CGI_COL));
% beijing_wcdma_lac_ci = cell(0);
% for ii=1:size(cgiList, 1)
%     GSMCellID = cgiList(ii);
%     lac = fix(GSMCellID/2^16);
%     cid = GSMCellID-lac*2^16;
%     
%     lac_ci = strcat(num2str(lac), '-', num2str(cid));
%     beijing_wcdma_lac_ci{ii, 1} = lac_ci;
% end

% cgiList = unique(measurement_Rural_Orangeville(:,CGI_COL));
% rural_orangeville_lac_ci = cell(0);
% for ii=1:size(cgiList, 1)
%     GSMCellID = cgiList(ii);
%     GSMCellID = GSMCellID+61106614173696;%+37937E080000
%     MCC = fix(GSMCellID/2^48)+284; %(2^56+2^52+2^51+2^50)/2^48
%     GSMCellID = GSMCellID - (MCC-284)*2^48;
%     GSMCellID = GSMCellID - fix(GSMCellID/2^44)*2^44;
%     MNC = fix(GSMCellID/2^32);
%     GSMCellID = GSMCellID-MNC * 2^32;
%     lac = fix(GSMCellID/2^16);
%     cid = GSMCellID-lac*2^16;
%     
%     lac_ci = strcat(num2str(lac), '-', num2str(cid));
%     rural_orangeville_lac_ci{ii, 1} = lac_ci;
% end

idx_lac_ci = 1;
idx_lat_long = 2;
idx_radius = 3;
lac_ci = cell_info{idx_lac_ci};

idx_start_angle = 1;
idx_stop_angle = 2;


if isempty(cell_angle)
    cell_angle = [0, 360];
end

start_angle = cell_angle(idx_start_angle);
stop_angle = cell_angle(idx_stop_angle);
[start_angle, stop_angle] = changeGeoAngle(start_angle, stop_angle);


%% rogers
% load rogers_msr;
% id_type = 2;
% measurementData = [];
% if find(ismember(urban_dt_lac_ci, lac_ci)>0)
%     measurementData = measurement_Urban_DT;
%     id_type = 1;
% elseif find(ismember(suburban_scarborough_lac_ci, lac_ci)>0)
%     measurementData = measurement_Suburban_Scarborough;
% elseif find(ismember(suburban_brampton_lac_ci, lac_ci)>0)
%     measurementData = measurement_Suburban_Brampton;
% elseif find(ismember(rural_orangeville_lac_ci, lac_ci)>0)
%     measurementData = measurement_Rural_Orangeville;
% end

%% turkey
% load turkey_msr;
% id_type = 1;
% measurementData = [];
% if find(ismember(beyoglu_2g_lac_ci, lac_ci)>0)
%     measurementData = measurement2G_beyoglu;
% elseif find(ismember(incek_2g_lac_ci, lac_ci)>0)
%     measurementData = measurement2G_incek;
% elseif find(ismember(sariyer_2g_lac_ci, lac_ci)>0)
%     measurementData = measurement2G_sariyer;
% elseif find(ismember(ulus_2g_lac_ci, lac_ci)>0)
%     measurementData = measurement2G_ulus;
% elseif find(ismember(uskudar_2g_lac_ci, lac_ci)>0)
%     measurementData = measurement2G_uskudar;
% end

%% beijing
load beijing_msr;
id_type = 3;
measurementData = measurement;

%% wuxi
% load wx_msr;
% id_type = 1;
% measurementData = measurement_134;

measurementData = sortrows(measurementData, [CGI_COL]);

ClusterSchemeLayer = 1;
tagMatrix = (1:1:length(measurementData(:,1)))'*ones(1,ClusterSchemeLayer);
for i = 1:length(measurementData(:,CGI_COL))-1
    % Layer 1
    if measurementData(i, CGI_COL) == measurementData(i+1, CGI_COL)
        tagMatrix(i+1, 1) = tagMatrix(i, 1);    
    end
end

cgiList = unique(measurementData(:,CGI_COL));
% numOfNB=MaxGsmNeighborCellCalculated;
for numOfCgi = 1:length(cgiList(:,1))
    for Layer = 1:ClusterSchemeLayer
        CGI = cgiList(numOfCgi,:);
        CgiStartPointIndex = find(CGI==measurementData(:,CGI_COL), 1, 'first'); %Extract points of same CGI
        CgiEndPointIndex = find(CGI==measurementData(:,CGI_COL), 1, 'last');
     
        GSMCellID = CGI;
        
        if id_type==2
            GSMCellID = GSMCellID+61106614173696;%+37937E080000
            MCC = fix(GSMCellID/2^48)+284; %(2^56+2^52+2^51+2^50)/2^48
            GSMCellID = GSMCellID - (MCC-284)*2^48;
            GSMCellID = GSMCellID - fix(GSMCellID/2^44)*2^44;
            MNC = fix(GSMCellID/2^32);
            GSMCellID = GSMCellID-MNC * 2^32;
        end
        % MCC=302;
        % MNC=720;
%         LAC = fix(GSMCellID/2^16);
%         CID = GSMCellID-LAC*2^16;
        
        lac = fix(GSMCellID/2^16);
        cid = GSMCellID-lac*2^16;
        
        
        lac_cid = strcat(num2str(lac), '-', num2str(cid)); 
        
        if id_type==3
            lac_cid = num2str(GSMCellID);
        end
        
        
        if strcmp(lac_cid, lac_ci)
            refPoint = deg2rad(cell_info{idx_lat_long});
        
%         sizeOfCluster = CgiEndPointIndex-CgiStartPointIndex+1;
            LatLongArrayOfCluster=measurementData(CgiStartPointIndex:CgiEndPointIndex,[LAT_COL LONG_COL INDEX_COL]);
            LatLongArrayOfCluster=sortrows(LatLongArrayOfCluster,3);
            LatLongArrayOfCluster=LatLongArrayOfCluster(:,1:2);

            LatLongArrayOfCluster = deg2rad(LatLongArrayOfCluster);

            enu_lat_long = convertlatlong2enu(LatLongArrayOfCluster, refPoint);
            
            % draw plot
            figure;
            hold on;
            plot(enu_lat_long(:, 1), enu_lat_long(:, 2), '*');
            
            % draw arc
            p.x = 0;
            p.y = 0;
            radius = cell_info{idx_radius};
            for jj=1:length(radius)
                drawArc(p, radius(jj), start_angle, stop_angle, 1);
                
                [X, Y ] = drawCircle(p.x, p.y, radius(jj));
                plot(X, Y, '--');
            end
            
            text(5, 5, lac_ci);
        end
    end
end

% end