% Script to make figures summarizing regression results
function main07_make_summary_figs(project,inference_id,use_weights,constrained_inference,ap_raw_flag)

% dataPaths
inputOutputPath = ['../../out/input_output/' project '/'];
inference_string = [inference_id '_wt' num2str(use_weights) '_apRaw' num2str(ap_raw_flag) '_ct' num2str(constrained_inference)];
regressionPath = [inputOutputPath '/results_' inference_string '/'];

% load data tables
load([inputOutputPath 'tf_input_struct.mat'])
regression_set = readtable([inputOutputPath 'final_regression_set.csv']);
load([regressionPath '/input_output_results_final.mat']);

% make fig path
figPath = ['../../fig/input_output/' project '/results_' inference_string '/'];
mkdir(figPath)

%%% Extract key regression stats from structure
pd_variables = reg_struct(1).reg_variables;
base_variables = reg_struct(1).base_variables;

% Build full regression set
ap_vec = tf_input_struct(1).ap_vec;
ap_field_eve2 = reg_struct(1).ap_field_eve2;
InterpGrid = tf_input_struct(1).InterpGrid;

% Motif info
motif_string_cell = reg_struct(1).motif_string_cell;
motif_sign_cell = reg_struct(1).motif_sign_cell;
motif_name_cell = reg_struct(1).motif_name_cell;

% Assign colors to motifs
cm = hsv(128);
inc = floor(128/(numel(motif_name_cell)));
motif_colors =  cm(1+(0:inc:inc*numel(motif_name_cell)),:);
n_predictors = numel(pd_variables);

%%% Build full regression set (more efficient to do this on the fly)
final_regression_set = add_tf_fields(regression_set,tf_input_struct,ap_field_eve2);


% update regression table to contain all necessary predictor fields
regression_inputs = augment_structure(final_regression_set, base_variables, pd_variables);
% define divergent colormap
cm = jet(128);
red = cm(120,:);
blue = cm(20,:);
rb_mat = [repmat(red,64,1) ; repmat(blue,64,1)];
white = repmat([1 1 1]/5*4,128,1);
shade_vec = [linspace(0,1,64)' ; flipud(linspace(0,1,64)')];
rb_cm = flipud(shade_vec.*white + (1-shade_vec).*rb_mat);

% %% Specify variables to examine
pt_vars = {'bin_act_flag','z_on_flag'};
lbl_strings = {'fraction active', 'fraction bursting'};
reg_flag_key = {reg_struct.dp_flag};


for i = 1:numel(pt_vars)    
    dp_var = reg_struct(i).dp_var;        
    dp_flag_var = reg_struct(i).dp_flag;
    % make directory
    VarPath = [figPath dp_flag_var '/'];
    mkdir(VarPath);
    
    use_weights = reg_struct(1).use_weights;
    dp_wt_vec = final_regression_set.([dp_flag_var '_wt']);
    if ~use_weights
        dp_wt_vec = ~isnan(dp_wt_vec);
    end    
    base_variables = reg_struct(i).base_variables;
    flag_vec = regression_inputs.(dp_flag_var);
    
    filt_set = regression_inputs(flag_vec==1&dp_wt_vec>0,:);
    
    filt_set.AP = round(filt_set.(ap_field_eve2)/10)*10;
    filt_set.Time = round(filt_set.Time/60)*60;
    % get averages TF inputs by AP and Time
    pd_array = grpstats(filt_set,{'AP','Time'},{'mean'},'DataVars',pd_variables);            
    % Extract results structure
    reg_results = reg_struct(i).reg_results;   
    
    % extract fit scores and n var data
    score_vec = [reg_results.score_test];    
    score_vec = score_vec - min(score_vec);
    n_var_vec = [reg_results.n_vars];
    motif_index_vector = [reg_results.motif_index];
    
    %%%%%%%%%%%%%% Make scatter plot of fit metric vs. N states %%%%%%%%%%%    
    % Color code dots according to motif    
    ft_scatter = figure('Visible','off');
    hold on    
    jitter = (rand(size(n_var_vec))-0.5)*.15;
    % No Motif assigned
    s = scatter(n_var_vec(motif_index_vector==0) + jitter(motif_index_vector==0)...
        ,score_vec(motif_index_vector==0),30,'MarkerFaceColor',[.3 .3 .3] ,...
                'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.5);    
    lgd_str = {'other'};
    for j = 1:numel(motif_string_cell)
        ft = motif_index_vector==j;
        if sum(ft) > 0
            s = [s scatter(n_var_vec(ft) + jitter(ft) ,score_vec(ft),30,'MarkerFaceColor', motif_colors(j,:),...
                    'MarkerEdgeAlpha',1,'MarkerFaceAlpha',1,'MarkerEdgeColor','black')];     
            lgd_str = [lgd_str{:} motif_name_cell(j)];
        end
    end
    xlabel('number of transcription factors')
    ylabel('relative log-likelihood (10^3)')    
    p = plot(0,0);    
    xlim([.75,4.25])
    ylim([0 1.05*max(score_vec)]);
    set(gca,'XTick',[1:4])
%     set(gca,'YTick',[0:2:10]*1E3)
%     set(gca,'YTickLabel',[0:2:10])
    set(gcf,'Name','Model Fit vs. Number of Predictors');    
%     lgd = legend(s,lgd_str{:},'Location','southoutside','Orientation','horizontal','NumColumns',3);
    lgd.FontSize = 10;%,'FontSize',10)    
    StandardFigurePBoC(p,gca)           
    saveas(ft_scatter, [VarPath 'score_scatter.png'])
        
    %%%%%%%%%%%%%%%% Make Coefficient Map Figure %%%%%%%%%%%%%%%%%%%%%%%%%    
    % Check to see how many times each variable is permitted to appear
    max_multi_vec = reg_struct(i).max_multi_vec;
    max_n = sum(max_multi_vec);
    % Generate string for x axis label
    x_lbl = {};
    for j = 1:numel(base_variables)
        for k = 1:max_multi_vec(j)
            x_lbl = [x_lbl{:} {[base_variables{j} num2str(k)]}];
        end
    end
    % initialize array to store coefficient info
    beta_array = NaN(numel(score_vec),max_n);   
    % obtain vector of vars by model, and corresponding coefficients. 
    [~, si] = sort(score_vec,'descend');
    ind = 1;   
    [~,sort_ref] = sort(pd_variables);
    for k = si                
        impact = reg_results(k).sign_vec.*reg_results(k).pd_contributions_norm;        
        impact = impact(2:end);        
        pdv = reg_results(k).pd_vars;
        for m = 1:numel(base_variables)
            matches = find(contains(pdv,base_variables{m}));
            indices = find(contains(x_lbl,base_variables{m}));
            for n = 1:numel(matches)
                beta_array(ind,indices(n)) = impact(matches(n));
            end
        end
        ind = ind+1;
    end            
    % figure summarizing model architectures
    beta_fig = figure('Visible','off');
    hold on
    colormap(rb_cm)
    beta_array_plot = flipud(beta_array);            

    p = imagesc(beta_array_plot);
    set(p,'AlphaData',~isnan(beta_array_plot))   
    caxis([-1 1])    
    set(gca,'xtick',1:n_predictors,'xticklabels',x_lbl);    
    xtickangle(-45)    
    set(gca,'ytick',[])
    h = colorbar;
    ylabel(h,'slope (\beta)')
    xlabel('transcription factor')
    ylabel('model')
    set(gcf,'Name','Predictor Coefficients Sorted By Model Fit');   
    axis([.5 size(beta_array_plot,2)+.5 0 size(beta_array_plot,1)])
    beta_fig.InvertHardcopy = 'off';
    beta_fig.Color = 'white';
    box on
    saveas(beta_fig, [VarPath 'coefficient_map.png'])    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Figures to Explore Impact of Predictor Power (Order) %%%%%%
    
    % Plot Average Parameter Order vs. Model Score...Are higher powers
    % favored in general?
    mean_power_vec = NaN(size(score_vec));
    for j = 1:numel(reg_results)
        pdv = reg_results(j).pd_vars;
        pwv = ones(size(pdv));
        for k = 1:numel(pdv)
            tf = pdv{k};
            pw = tf(strfind(tf,'_')+1:end);
            if ~isempty(pw)
                pwv(k) = str2double(pw);
            end
        end
        mean_power_vec(j) = mean(pwv);
    end
    
    power_scatter = figure('Visible','off');
    hold on        
    % No Motif assigned
    jitter = (rand(size(n_var_vec))-0.5)*.1;
    ft = motif_index_vector==0;
    s = scatter(mean_power_vec(ft) + jitter(ft)...
        ,score_vec(ft),30,'MarkerFaceColor',[.3 .3 .3] ,...
                'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);    
    % iterate through motifs
    lgd_str = {'other'};
    for j = 1:numel(motif_string_cell)
        ft = motif_index_vector==j;
        if sum(ft) > 0
            s = [s scatter(mean_power_vec(ft)+jitter(ft),score_vec(ft),30,'MarkerFaceColor', motif_colors(j,:),...
                    'MarkerEdgeAlpha',1,'MarkerFaceAlpha',.7,'MarkerEdgeColor','black')];     
            lgd_str = [lgd_str{:} motif_name_cell(j)];
        end
    end
    xlabel('number of transcription factors')
    ylabel('relative log-likelihood (10^3)')    
    p = plot(0,0);
%     lgd = legend(s,lgd_str{:},'Location','southoutside','Orientation','horizontal','NumColumns',3);
    lgd.FontSize = 10;%,'FontSize',10)     
    xlim([.75,1.1*max(mean_power_vec)])
    ylim([0 1.05*max(score_vec)]);
%     set(gca,'XTick',[1:4])
%     set(gca,'YTick',[0:2:10]*1E3)
%     set(gca,'YTickLabel',[0:2:10])
    set(gcf,'Name','Model Fit vs. Mean Predictor Order'); 
    grid on
%     StandardFigurePBoC(p,gca)        
    saveas(power_scatter, [VarPath 'power_score_scatter.png'])
    
    % Now loop through individual factors to see there is a preference for
    % power on the single variable level. Highloght models for which var of
    % interest is actually relevant
    sym_colors = {blue [.5 .5 .5] red};
    sym_names = {'other (neg)', 'other (bi)','other (pos)'};
    for j = 1:numel(base_variables)
        pd_var = base_variables{j};
        % iterate through each model and identify var power, signficicance,
        % and function (activating or repressing). If variable is not
        % present in model, power assigned as 0
        mean_var_power = zeros(1,numel(reg_results));
        var_effect = NaN(1,numel(reg_results));
        for k = 1:numel(reg_results)
            matches = contains(reg_results(k).pd_vars,pd_var);
            wt_vec = reg_results(k).pd_contributions_norm;
            wt_vec = wt_vec(2:end);
            wt_vec = wt_vec(matches);
            % exrract powers
            pdv = reg_results(k).pd_vars(matches);
            pwv = ones(size(pdv));
            for m = 1:numel(pdv)
                tf = pdv{m};
                pw = tf(strfind(tf,'_')+1:end);
                if ~isempty(pw)
                    pwv(m) = str2double(pw);
                end
            end
            if sum(reg_results(k).sig_vec(matches)) > 0                                
                mean_var_power(k) = sum(pwv.*wt_vec) / sum(wt_vec);                
                sign_vec = reg_results(k).sign_vec(2:end);
                var_effect(k) = sum(sign_vec(matches).*reg_results(k).sig_vec(matches));                
            % if no effects are significant, just pull and average powers
            elseif sum(matches) > 0
                mean_var_power(k) = mean(pwv);
            end
        end
        var_sig_vec = ~isnan(var_effect);
        % only display var effect when no corresponding motif exists
        var_effect(motif_index_vector~=0) = NaN;
        
        % now plot 
        jitter = (rand(size(n_var_vec))-0.5)*.1;
        var_power_fig = figure('Visible','off');
        hold on        
        % No Motif assigned
        ft = motif_index_vector==0;
        s = scatter(mean_var_power(ft)+jitter(ft)...
            ,score_vec(ft),30,'MarkerFaceColor',[.3 .3 .3] ,...
                    'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.3);    
        % iterate through motifs
        lgd_str = {'other'};
        for k = 1:numel(motif_string_cell)
            ft = motif_index_vector==k;
            ms = 20;
            ma = .3;
            ea = 0;
            if max(var_sig_vec(ft))==1
                ms = 30;
                ma = .7;
                ea = 1;
            end
            if sum(ft) > 0
                s = [s scatter(mean_var_power(ft)+jitter(ft),score_vec(ft),ms,'MarkerFaceColor', motif_colors(k,:),...
                        'MarkerEdgeAlpha',ea,'MarkerFaceAlpha',ma,'MarkerEdgeColor','black')];     
                lgd_str = [lgd_str{:} motif_name_cell(k)];
            end
        end
        % lastly, plot models where variable had an impact and which did
        % not fall into a motif bin
        eff_index = unique(var_effect(~isnan(var_effect)));
        for k = 1:numel(eff_index)
            ft = var_effect == eff_index(k);
            s = [s scatter(mean_var_power(ft)+jitter(ft),score_vec(ft),30,'^','MarkerFaceColor', sym_colors{eff_index(k)+2},...
                        'MarkerEdgeAlpha',1,'MarkerFaceAlpha',1,'MarkerEdgeColor','black')];     
            lgd_str = [lgd_str{:} sym_names(eff_index(k)+2)];  
        end
        
        xlabel(['predictor order (' pd_var ')'])
        ylabel('relative log-likelihood (10^3)')    
        p = plot(0,0);
%         lgd = legend(s,lgd_str{:},'Location','southoutside','Orientation','horizontal','NumColumns',3);
        lgd.FontSize = 10;%,'FontSize',10)   
        xlim([-.25, 1.1*max(mean_var_power)])
        ylim([0 1.25*max(score_vec)]);
%         set(gca,'XTick',[0:5])
%         set(gca,'YTick',[0:2:10]*1E3)
%         set(gca,'YTickLabel',[0:2:10])
        set(gcf,'Name',['Model Fit vs. Mean Predictor Order (' pd_var ')']);    
%         StandardFigurePBoC(p,gca)                   
        saveas(var_power_fig, [VarPath 'power_scatter_' pd_var '.png']);
    end    
    %%%%%%%%%%% Decompose Model Error by Time and AP %%%%%%%%%%%%%%%%%%%%%%
    ap_lims = [min(filt_set.AP) max(filt_set.AP)];
    ap_cal_factor = pi/(ap_lims(2)-ap_lims(1))/2;
    time_lims = [min(filt_set.Time) max(filt_set.Time)];
    time_cal_factor = pi/(time_lims(2)-time_lims(1))/2;
    % extract error vectors
    ap_err_vec = [reg_results.ap_err];
    ap_angle_vec = pi/2 - (ap_err_vec - ap_lims(1))*ap_cal_factor;
    time_err_vec = [reg_results.time_err];
    time_angle_vec = pi/2 - (time_err_vec - time_lims(1))*time_cal_factor;
    
    % Make AP scatter
    ap_score_fig = figure('Visible','off');
    a_score_vec = score_vec.*cos(ap_angle_vec);
    p_score_vec = score_vec.*sin(ap_angle_vec);
    hold on
    % No Motif assigned
    s = scatter(a_score_vec(motif_index_vector==0)...
        ,p_score_vec(motif_index_vector==0),30,'MarkerFaceColor',[.3 .3 .3] ,...
                'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.5);    
    lgd_str = {'other'};
    for j = 1:numel(motif_string_cell)
        ft = motif_index_vector==j;
        if sum(ft) > 0
            s = [s scatter(a_score_vec(ft) ,p_score_vec(ft),30,'MarkerFaceColor', motif_colors(j,:),...
                    'MarkerEdgeAlpha',1,'MarkerFaceAlpha',1,'MarkerEdgeColor','black')];     
            lgd_str = [lgd_str{:} motif_name_cell(j)];
        end
    end
    xlabel('anterior score')
    ylabel('posterior score')    
    p = plot(0,0);
%     lgd = legend(s,lgd_str{:},'Location','southoutside','Orientation','horizontal','NumColumns',3);
    lgd.FontSize = 10;%,'FontSize',10)    
    xlim([0 1.05*max([a_score_vec p_score_vec])])
    ylim([0 1.05*max([a_score_vec p_score_vec])]);
%     set(gca,'XTick',[0:2:10]*1E3)
%     set(gca,'YTick',[0:2:10]*1E3)
%     set(gca,'YTickLabel',[0:2:10])
    set(gcf,'Name','Model Fit vs. AP');    
%     StandardFigurePBoC(p,gca)   
    saveas(ap_score_fig, [VarPath 'ap_score_scatter.png'])
    
    % Make Time scatter
    time_score_fig = figure('Visible','off');
    e_score_vec = score_vec.*cos(time_angle_vec);
    l_score_vec = score_vec.*sin(time_angle_vec);
    hold on
    % No Motif assigned
    s = scatter(e_score_vec(motif_index_vector==0)...
        ,l_score_vec(motif_index_vector==0),30,'MarkerFaceColor',[.3 .3 .3] ,...
                'MarkerEdgeAlpha',0,'MarkerFaceAlpha',.5);    
    lgd_str = {'other'};
    for j = 1:numel(motif_string_cell)
        ft = motif_index_vector==j;
        if sum(ft) > 0
            s = [s scatter(e_score_vec(ft) ,l_score_vec(ft),30,'MarkerFaceColor', motif_colors(j,:),...
                    'MarkerEdgeAlpha',1,'MarkerFaceAlpha',1,'MarkerEdgeColor','black')];     
            lgd_str = [lgd_str{:} motif_name_cell(j)];
        end
    end
    xlabel('early score')
    ylabel('late score')    
%     p = plot(0,0);
%     lgd = legend(s,lgd_str{:},'Location','southoutside','Orientation','horizontal','NumColumns',3);
    lgd.FontSize = 10;%,'FontSize',10)    
    xlim([0 1.05*max([a_score_vec p_score_vec])])
    ylim([0 1.05*max([a_score_vec p_score_vec])]);
%     set(gca,'XTick',[0:2:10]*1E3)
%     set(gca,'YTick',[0:2:10]*1E3)
%     set(gca,'YTickLabel',[0:2:10])
    set(gcf,'Name','Model Fit vs. Time');    
%     StandardFigurePBoC(p,gca)   
    saveas(time_score_fig, [VarPath 'time_score_scatter.png'])    
    
    %%%%%%%%%% Heat Maps comparing model predictions to actual response%%%%
    % Generate maps for all motifs, irrespective of score
    % Generate maps for max(top 10%,20) of models for each N 
    
    % Make Heatmap depicting true response
    ap_grp = pd_array.AP;
    time_grp = pd_array.Time;
    ap_all = filt_set.AP;
    time_all = filt_set.Time;
            
    ap_ref = unique(ap_grp);
    time_ref = unique(time_grp);
    ap_num = numel(ap_ref);
    t_num = numel(time_ref);
    
    y_grid = NaN(numel(time_ref),numel(ap_ref));
    dp_var_vec = filt_set.(dp_var);
    for n = 1:numel(ap_grp)
        y_grid(time_grp(n)==time_ref,ap_grp(n)==ap_ref) = nanmean(dp_var_vec(ap_all==ap_grp(n)&time_all==time_grp(n)));
    end
    % Make figure
    ub = prctile(y_grid(:),99);
    lb = prctile(y_grid(:),1);
    
    true_fig = figure('Visible','off');
    colormap(jet(128));
    % actual    
    imagesc(y_grid);            
    set(gca,'xtick',1:2:ap_num,'xticklabels',round(ap_ref(1:2:ap_num)/10));
    set(gca,'ytick',1:2:t_num,'yticklabels',round(time_ref(1:2:t_num)/60));                       
    xtickangle(-45)
    h = colorbar;
    caxis([lb ub])
    set(gcf,'Name',['Actual Response: ' dp_flag_var]);      
    ylabel('minutes')
    ylabel(h, lbl_strings{i})
    xlabel('AP')
    saveas(true_fig,[VarPath 'response_heatmap.png'])
    
    % Make Heatmaps for qualifying models
    pd_col_ref =  pd_array.Properties.VariableNames;
    for n = 1:numel(pd_col_ref)
        pd_col_ref{n} = pd_col_ref{n}(strfind(pd_col_ref{n},'_')+1:end);
    end    
    hm_vec = false(size(motif_index_vector));
    for n = 1:max(n_var_vec)
        n_indices = find(n_var_vec==n);
        p10 = prctile(score_vec(n_indices),90);
        ft = n_var_vec==n&score_vec>=p10;
        if numel(n_indices) < 20
            hm_vec = hm_vec | n_var_vec==n;  
        else
            score_temp = score_vec;
            score_temp(n_var_vec~=n) = -Inf;            
            [~, si] = sort(score_temp,'descend');
            hm_vec(si(1:20)) = true;
        end
    end
    
    % Generate Prediction heatmaps
    pd_hm_array = NaN(size(y_grid,1),size(y_grid,2),sum(hm_vec));
    cf_hm_array = NaN(size(y_grid,1),size(y_grid,2),sum(hm_vec));
    iter = 1;
    mkdir([VarPath '/pd_profiles/'])
    hm_indices = find(hm_vec);
    % Generate heatmaps for all results    
    for k = hm_indices
        pd_vars = reg_results(k).pd_vars;
        var_array = [];
        for n = 1:numel(pd_vars)
            var_array = [var_array pd_array{:,ismember(pd_col_ref,pd_vars{n})}];
        end        
        beta_fit = reg_results(k).beta_fit;
        pd_vec = mnrval(beta_fit,var_array);   
        pd_grid = NaN(size(y_grid));           
        for n = 1:numel(ap_grp)            
            pd_grid(time_grp(n)==time_ref,ap_grp(n)==ap_ref) = pd_vec(n,2);
        end        
        cf_grid = pd_grid - y_grid;
        pd_hm_array(:,:,iter) = pd_grid;
        cf_hm_array(:,:,iter) = cf_grid;
        iter = iter + 1;
        
        % make save string
        s_string = pd_vars{1}; 
        for n = 2:numel(pd_vars)
            s_string = [s_string '_' pd_vars{n}];
        end
        % Make figure
        ub = prctile(y_grid(:),99);
        lb = prctile(y_grid(:),1);
        
        pd_fig = figure('Visible','off');
        colormap(jet(128));
        % actual    
        imagesc(pd_grid);            
        set(gca,'xtick',1:2:ap_num,'xticklabels',round(ap_ref(1:2:ap_num)/10));
        set(gca,'ytick',1:2:t_num,'yticklabels',round(time_ref(1:2:t_num)/60));                       
        xtickangle(-45)
        h = colorbar;
%         caxis([lb ub])
        if motif_index_vector(k) ~= 0
            title(motif_name_cell{motif_index_vector(k)});
        else
            title('other')
        end
%         title(['Predicted Response: ' lbl_strings{j} ' (' s_string ')'],'interpreter','none');      
        ylabel('minutes')
        ylabel(h, lbl_strings{i})
        xlabel('AP')

        saveas(pd_fig, [VarPath '/pd_profiles/N' num2str(n_var_vec(k)) ...
            '_S' num2str(round(score_vec(k))) ...
            '_' s_string '.png'])
        saveas(pd_fig, [VarPath '/pd_profiles/N' num2str(n_var_vec(k)) ...
            '_S' num2str(round(score_vec(k))) ...
            '_' s_string '.fig'])
        close all
    end
    
    % Generate confusion heatmaps    
    iter = 1;
    mkdir([VarPath '/cf_profiles/'])
    hm_indices = find(hm_vec);    
    % Generate heatmaps for all results
    for k = hm_indices                
        pd_vars = reg_results(k).pd_vars;
        cf_grid = cf_hm_array(:,:,iter);        
        % make save string
        s_string = pd_vars{1}; 
        for n = 2:numel(pd_vars)
            s_string = [s_string '_' pd_vars{n}];
        end
        % Make figure
        ub = 1;
        lb = -1;        
        iter = iter + 1;
        pd_fig = figure('Visible','off');
        colormap(rb_cm);
        % actual    
        imagesc(cf_grid);     
        if motif_index_vector(k) ~= 0
            title(motif_name_cell{motif_index_vector(k)});
        else
            title('other')
        end
        set(gca,'xtick',1:2:ap_num,'xticklabels',round(ap_ref(1:2:ap_num)/10));
        set(gca,'ytick',1:2:t_num,'yticklabels',round(time_ref(1:2:t_num)/60));                       
        xtickangle(-45)
        h = colorbar;
%         caxis([lb ub])
%         title(['Predicted Response: ' lbl_strings{j} ' (' s_string ')'],'interpreter','none');      
        ylabel('minutes')
        ylabel(h, lbl_strings{i})
        xlabel('AP')

        saveas(pd_fig, [VarPath '/cf_profiles/N' num2str(n_var_vec(k)) ...
            '_S' num2str(round(score_vec(k))) ...
            '_' s_string '.png'])          
    end
end



