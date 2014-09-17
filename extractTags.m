
% extract the tags

% beijing tags
% for ii=1:length(beijing_tags)
%     tag = beijing_tags{ii, 2};
%     
%     if isempty(strfind(tag, ','))
%         % just cell tag
%         beijing_cell_tags= [beijing_cell_tags; beijing_tags(ii, :)];
%     end
% end

% beijing_wcdma_cells = [];
% beijing_wcdma_tags = [];
% [~, idx] = unique(Beijing_cell(:, 14));
% unq_beijing_cell = Beijing_cell(idx, :);
% for ii=1:length(unq_beijing_cell)
%     cell_id = unq_beijing_cell(ii, 14);
%     lac = fix(cell_id/2^16);
%     cid = mod(cell_id, 2^16);
%     
%     lac_ci = strcat(num2str(lac), '-', num2str(cid)); 
%     
%     if ~isempty(beijing_wcdma_cells) && ~isempty(find(beijing_wcdma_cells(:, 14)==cell_id))
%         continue;
%     end
%     
%     for jj=1:length(beijing_cell_tags)
%         cell2_string = beijing_cell_tags{jj, 2};
%         cell2_split = regexp(cell2_string, '-', 'split');
%         
%         lac_ci2 = strcat(cell2_split{3}, '-', cell2_split{4});
%         
%         if lac_ci==lac_ci2
%             beijing_wcdma_cells = [beijing_wcdma_cells; Beijing_cell(ii, :)];
%             beijing_wcdma_tags = [beijing_wcdma_tags; beijing_cell_tags(jj, :)];
%             break;
%         end
%     end
% end

% wuxi tags

% wx_tags = [wx_tags; wx_134_tags];
% wx_wcdma_cells = [];
% wx_wcdma_tags = [];
% 
% for ii=1:length(wx_cell)
%     cid_1 = wx_cell(ii, 1);
%     
% %     if ~isempty(wx_wcdma_cells) && ~isempty(find(wx_wcdma_cells(:, 1)==cid_1))
% %         continue;
% %     end
%     
%     for jj=1:length(wx_tags)
%         cell2_string = wx_tags{jj, 2};
%         cell2_split = regexp(cell2_string, '-', 'split');
%         
%         cid_2 = str2num(cell2_split{4});
%         
%         if cid_1==cid_2
%             wx_wcdma_cells = [wx_wcdma_cells; wx_cell(ii, :)];
%             wx_wcdma_tags = [wx_wcdma_tags; wx_tags(jj, :)];
%             break;
%         end
%     end
% end

% remove NaN cell for wuxi
wx_new_cell = [];
for ii=1:length(wx_cell)
    if any(isnan(wx_cell(ii, :)))
        continue;
    end
    
    wx_new_cell = [wx_new_cell; wx_cell(ii, :)];
end