function [] = function_iden_strategy(SA_flag,totkk,Corr_matrix,threshold,gsax_labels)
if totkk~=1
    add_Corr=load(['SA_re_',num2str(totkk-1),'.mat']);
    add_Corr=add_Corr.iden_par;   
    temp1=fliplr(gsax_labels);
    index_row=[];
    for jj=1:1:size(add_Corr.tot_row,2)
        value_to_find=add_Corr.tot_row{jj};
       indices = find(cellfun(@(x) isequal(x, value_to_find), temp1));
       index_row=[index_row,indices];
    end       
    temp2=(gsax_labels);
    index_col=[];
    for jj=1:1:size(add_Corr.tot_col,2)
        value_to_find=add_Corr.tot_col{jj};
       indices = find(cellfun(@(x) isequal(x, value_to_find), temp2));
       index_col=[index_col,indices];
    end   
    Corr_matrix(:,index_row)=-0.2;
    Corr_matrix(index_col,:)=-0.2;    
end

Iden_par_SA=[];
corr_index=[];
for kk=13:-1:1
    corr_ref=Corr_matrix(:,kk);
    if kk==13
        corr_index_temp1=find(corr_ref>=threshold);
        corr_index_temp2=13-find(corr_ref>=threshold)+1;
        corr_index=union(corr_index_temp2,corr_index);
        Corr_matrix(corr_index_temp1,:)=NaN;
        Corr_matrix(:,corr_index_temp2)=NaN;
    else
        corr_ref=Corr_matrix(:,kk);      
        corr_index_temp1=find(corr_ref>=threshold);
        corr_index_temp2=13-find(corr_ref>=threshold)+1;        
        Corr_matrix(corr_index_temp1,:)=NaN;
        Corr_matrix(:,corr_index_temp2)=NaN;
    end   
end

% %% Create figures for each identification step
% % 1. Model parameters identified in this step are shown in yellow.
% % 2. Model parameters obtained from the last step are shown in purple.
% % 3. Model parameters that are correlated with each other are shown in black.
% 
% % Create an figure window
% hFig = figure('Visible', 'on');
% 
% % Create the heatmap
% hHeatmap = heatmap(Corr_matrix);
% hHeatmap.XDisplayLabels = fliplr(gsax_labels);
% hHeatmap.YDisplayLabels = gsax_labels;
% 
% % Get and set the colormap
% cmap = colormap(hHeatmap);
% yellow_color = [1, 1, 0]; 
% purple_color = [0.5, 0, 0.5]; 
% color_limits = hHeatmap.ColorLimits;
% num_colors = size(cmap, 1);
% yellow_indices = find(Corr_matrix > -0.1 & Corr_matrix <= 0);
% purple_indices = find(Corr_matrix <= -0.1);
% 
% % Set yellow regions
% for idx = yellow_indices'
%     [row, col] = ind2sub(size(Corr_matrix), idx);
%     cmap_index = round((Corr_matrix(row, col) - color_limits(1)) / (color_limits(2) - color_limits(1)) * (num_colors - 1)) + 1;
%     cmap(cmap_index, :) = yellow_color;
% end
% 
% % Set purple regions
% for idx = purple_indices'
%     [row, col] = ind2sub(size(Corr_matrix), idx);
%     cmap_index = round((Corr_matrix(row, col) - color_limits(1)) / (color_limits(2) - color_limits(1)) * (num_colors - 1)) + 1;
%     cmap(cmap_index, :) = purple_color;
% end
% 
% % Apply the customized colormap to the heatmap
% colormap(hHeatmap, cmap);
% caxis(hHeatmap, color_limits);
% xlabel(hHeatmap, 'ESPM parameters');
% ylabel(hHeatmap, 'ESPM parameters');
% hHeatmap.FontSize = 14;
% hHeatmap.CellLabelColor = 'none';
% set(hFig, 'unit', 'centimeters', 'position', [13 5 18 15]);
% 
% % Save the figure as a .fig file
% saveas(hFig, ['Iden_strategy',num2str(totkk),'.fig']);
% 
% % Close the figure window
% close(hFig);
% 
% 
% [raw_row,raw_col]=find(Corr_matrix==0);
% temp1=fliplr(gsax_labels);
% temp2=gsax_labels;
% iden_par.row=temp2(raw_row);
% iden_par.col=temp1(raw_col);
% if totkk==1
% iden_par.tot_row=iden_par.row;
% iden_par.tot_col=iden_par.col;
% elseif totkk~=1
% pre_iden_par=load(['SA_re_',num2str(totkk-1),'.mat']);  
% pre_iden_par=pre_iden_par.iden_par;
% iden_par.tot_row=union(iden_par.row,pre_iden_par.row);
% iden_par.tot_col=union(iden_par.col,pre_iden_par.col);
% end


% -------------------------------------------------------------------------
[raw_row,raw_col]=find(Corr_matrix==0);
temp1=fliplr(gsax_labels);
temp2=gsax_labels;
iden_par.row=temp2(raw_row);
iden_par.col=temp1(raw_col);
if totkk==1
iden_par.tot_row=iden_par.row;
iden_par.tot_col=iden_par.col;
elseif totkk~=1
pre_iden_par=load(['SA_re_',num2str(totkk-1),'.mat']);  
pre_iden_par=pre_iden_par.iden_par;
iden_par.tot_row=union(iden_par.row,pre_iden_par.row);
iden_par.tot_col=union(iden_par.col,pre_iden_par.col);
end
% -------------------------------------------------------------------------



save(['SA_re_',num2str(SA_flag),'.mat'],'iden_par')
end

