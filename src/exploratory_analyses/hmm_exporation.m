clear
close all
dataPath = 'E:\Annika\mse_dissection\out\ss_inference_results\mse_comparison_lateralML\';
load([dataPath 'hmm_results_summary_Gt'])
gt_hmm_results = hmm_results;
load([dataPath 'hmm_results_summary_Wt'])
wt_hmm_results = hmm_results;
load([dataPath 'hmm_results_summary_Hb'])
hb_hmm_results = hmm_results;
%% effective on rate comparisons
on_rate_fig = figure;
errorbar([wt_hmm_results.binID], [wt_hmm_results.eff_on_rate], [wt_hmm_results.eff_on_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([gt_hmm_results.binID], [gt_hmm_results.eff_on_rate], [gt_hmm_results.eff_on_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([hb_hmm_results.binID], [hb_hmm_results.eff_on_rate], [hb_hmm_results.eff_on_ste], 'Marker', '.', 'LineWidth', 1.25)
hold off

legend('Wt','Gt','Hb')
xlabel('AP Position (% Embryo Length)')
ylabel('Effective On Rate')
title('Effective On Rate Comparison')

saveas(on_rate_fig,'E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\eff_on_fig.png')
%% effective off rate comparisons
off_rate_fig = figure;
errorbar([wt_hmm_results.binID], [wt_hmm_results.eff_off_rate], [wt_hmm_results.eff_off_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([gt_hmm_results.binID], [gt_hmm_results.eff_off_rate], [gt_hmm_results.eff_off_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([hb_hmm_results.binID], [hb_hmm_results.eff_off_rate], [hb_hmm_results.eff_off_ste], 'Marker', '.', 'LineWidth', 1.25)
hold off

legend('Wt','Gt','Hb')
xlabel('AP Position (% Embryo Length)')
ylabel('Effective Off Rate')
title('Effective Off Rate Comparison')

saveas(off_rate_fig,'E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\eff_off_fig.png')

%% effective initiation rate comparisons
init_rate_fig = figure;
errorbar([wt_hmm_results.binID], [wt_hmm_results.eff_init_rate], [wt_hmm_results.eff_init_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([gt_hmm_results.binID], [gt_hmm_results.eff_init_rate], [gt_hmm_results.eff_init_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([hb_hmm_results.binID], [hb_hmm_results.eff_init_rate], [hb_hmm_results.eff_init_ste], 'Marker', '.', 'LineWidth', 1.25)
hold off

legend('Wt','Gt','Hb')
xlabel('AP Position (% Embryo Length)')
ylabel('Effective Initiation Rate')
title('Effective Initiation Rate Comparison')

saveas(init_rate_fig,'E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\eff_init_fig.png')

%% mean initiation ratio
init_ratio_fig = figure;
errorbar([wt_hmm_results.binID], [wt_hmm_results.init_ratio_mean], [wt_hmm_results.init_ratio_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([gt_hmm_results.binID], [gt_hmm_results.init_ratio_mean], [gt_hmm_results.init_ratio_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([hb_hmm_results.binID], [hb_hmm_results.init_ratio_mean], [hb_hmm_results.init_ratio_ste], 'Marker', '.', 'LineWidth', 1.25)
hold off

legend('Wt','Gt','Hb')
xlabel('AP Position (% Embryo Length)')
ylabel('Mean Initiation Ratio')
title('Mean Initiation Rate Comparison')

saveas(init_ratio_fig,'E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\mean_init_ratio_fig.png')

%% on ratio mean
on_ratio_fig = figure;
errorbar([wt_hmm_results.binID], [wt_hmm_results.on_ratio_mean], [wt_hmm_results.on_ratio_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([gt_hmm_results.binID], [gt_hmm_results.on_ratio_mean], [gt_hmm_results.on_ratio_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([hb_hmm_results.binID], [hb_hmm_results.on_ratio_mean], [hb_hmm_results.on_ratio_ste], 'Marker', '.', 'LineWidth', 1.25)
hold off

legend('Wt','Gt','Hb')
xlabel('AP Position (% Embryo Length)')
ylabel('Mean On Ratio')
title('Mean On Ratio Comparison')

saveas(on_ratio_fig,'E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\on_ratio_fig.png')

%% off ratio mean
off_ratio_fig = figure;
errorbar([wt_hmm_results.binID], [wt_hmm_results.off_ratio_mean], [wt_hmm_results.off_ratio_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([gt_hmm_results.binID], [gt_hmm_results.off_ratio_mean], [gt_hmm_results.off_ratio_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([hb_hmm_results.binID], [hb_hmm_results.off_ratio_mean], [hb_hmm_results.off_ratio_ste], 'Marker', '.', 'LineWidth', 1.25)
hold off

legend('Wt','Gt','Hb')
xlabel('AP Position (% Embryo Length)')
ylabel('Mean Off Ratio')
title('Mean Off Ratio Comparison')

saveas(off_ratio_fig,'E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\off_ratio_fig.png')

%% on occupancy effective mean
on_occ_eff_mean_fig = figure;
errorbar([wt_hmm_results.binID], [wt_hmm_results.on_occ_eff_mean], [wt_hmm_results.on_occ_eff_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([gt_hmm_results.binID], [gt_hmm_results.on_occ_eff_mean], [gt_hmm_results.on_occ_eff_ste], 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([hb_hmm_results.binID], [hb_hmm_results.on_occ_eff_mean], [hb_hmm_results.on_occ_eff_ste], 'Marker', '.', 'LineWidth', 1.25)
hold off

legend('Wt','Gt','Hb')
xlabel('AP Position (% Embryo Length)')
ylabel('Mean Effective Occupancy')
title('Mean Effective Occupancy Comparison')

saveas(on_occ_eff_mean_fig,'E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\on_occ_eff_mean_fig.png')

%% initiation mean
wt_initiation_mean = [wt_hmm_results.initiation_mean];
gt_initiation_mean = [gt_hmm_results.initiation_mean];
hb_initiation_mean = [hb_hmm_results.initiation_mean];
wt_initiation_std = [wt_hmm_results.initiation_std];
gt_initiation_std = [gt_hmm_results.initiation_std];
hb_initiation_std = [hb_hmm_results.initiation_std];
for i = 1:3
    initiation_mean_fig = figure;
    errorbar([wt_hmm_results.binID], wt_initiation_mean(i,:), wt_initiation_std(i,:), 'Marker', '.', 'LineWidth', 1.25)
    hold on
    errorbar([gt_hmm_results.binID], gt_initiation_mean(i,:), gt_initiation_std(i,:), 'Marker', '.', 'LineWidth', 1.25)
    hold on
    errorbar([hb_hmm_results.binID], hb_initiation_mean(i,:), hb_initiation_std(i,:), 'Marker', '.', 'LineWidth', 1.25)
    hold off

    legend('Wt','Gt','Hb')
    xlabel('AP Position (% Embryo Length)')
    ylabel('Mean Initiation')
    title(['Mean Initiation ' num2str(i) ' Comparison'])

    saveas(initiation_mean_fig,['E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\initiation_mean_fig' num2str(i) '.png'])
    
end

%% occupancy mean

wt_occupancy_mean = [wt_hmm_results.occupancy_mean];
gt_occupancy_mean = [gt_hmm_results.occupancy_mean];
hb_occupancy_mean = [hb_hmm_results.occupancy_mean];
wt_occupancy_std = [wt_hmm_results.occupancy_std];
gt_occupancy_std = [gt_hmm_results.occupancy_std];
hb_occupancy_std = [hb_hmm_results.occupancy_std];
for i = 1:3
    occupancy_mean_fig = figure;
    errorbar([wt_hmm_results.binID], wt_occupancy_mean(i,:), wt_occupancy_std(i,:), 'Marker', '.', 'LineWidth', 1.25)
    hold on
    errorbar([gt_hmm_results.binID], gt_occupancy_mean(i,:), gt_occupancy_std(i,:), 'Marker', '.', 'LineWidth', 1.25)
    hold on
    errorbar([hb_hmm_results.binID], hb_occupancy_mean(i,:), hb_occupancy_std(i,:), 'Marker', '.', 'LineWidth', 1.25)
    hold off

    legend('Wt','Gt','Hb')
    xlabel('AP Position (% Embryo Length)')
    ylabel('Mean Occupancy')
    title(['Mean Occupancy ' num2str(i) ' Comparison'])

    saveas(occupancy_mean_fig,['E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\occupancy_mean_fig' num2str(i) '.png'])
    
end

%% dwell mean

wt_dwell_mean = [wt_hmm_results.dwell_mean];
gt_dwell_mean = [gt_hmm_results.dwell_mean];
hb_dwell_mean = [hb_hmm_results.dwell_mean];
wt_dwell_std = [wt_hmm_results.dwell_std];
gt_dwell_std = [gt_hmm_results.dwell_std];
hb_dwell_std = [hb_hmm_results.dwell_std];
for i = 1:3
    dwell_mean_fig = figure;
    errorbar([wt_hmm_results.binID], wt_dwell_mean(i,:), wt_dwell_std(i,:), 'Marker', '.', 'LineWidth', 1.25)
    hold on
    errorbar([gt_hmm_results.binID], gt_dwell_mean(i,:), gt_dwell_std(i,:), 'Marker', '.', 'LineWidth', 1.25)
    hold on
    errorbar([hb_hmm_results.binID], hb_dwell_mean(i,:), hb_dwell_std(i,:), 'Marker', '.', 'LineWidth', 1.25)
    hold off

    legend('Wt','Gt','Hb')
    xlabel('AP Position (% Embryo Length)')
    ylabel('Mean Dwell')
    title(['Mean Dwell ' num2str(i) ' Comparison'])

    saveas(dwell_mean_fig,['E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\dwell_mean_fig' num2str(i) '.png'])
    
end

close all
%%
%plot second element of the R_fit_mean (transition from off state to middle
%state)

wt_rfit_mean = nan(1,numel(wt_hmm_results)); 
gt_rfit_mean = nan(1,numel(gt_hmm_results));
hb_rfit_mean = nan(1,numel(hb_hmm_results));

wt_rfit_std = nan(1,numel(wt_hmm_results)); 
gt_rfit_std = nan(1,numel(gt_hmm_results));
hb_rfit_std = nan(1,numel(hb_hmm_results));

for i = 1:numel(wt_rfit_mean)
    wt_rfit_mean(i) = wt_hmm_results(i).R_fit_mean(2);
    wt_rfit_std(i) = wt_hmm_results(i).R_fit_std(2);
end
for i = 1:numel(gt_rfit_mean)
    gt_rfit_mean(i) = gt_hmm_results(i).R_fit_mean(2);
    gt_rfit_std(i) = gt_hmm_results(i).R_fit_std(2);
end
for i = 1:numel(hb_rfit_mean)
    hb_rfit_mean(i) = hb_hmm_results(i).R_fit_mean(2);
    hb_rfit_std(i) = hb_hmm_results(i).R_fit_std(2);
end

rfit_fig = figure;
errorbar([wt_hmm_results.binID], wt_rfit_mean, wt_rfit_std, 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([gt_hmm_results.binID], gt_rfit_mean, gt_rfit_std, 'Marker', '.', 'LineWidth', 1.25)
hold on
errorbar([hb_hmm_results.binID], hb_rfit_mean, hb_rfit_std, 'Marker', '.', 'LineWidth', 1.25)
hold off

legend('Wt','Gt','Hb')
xlabel('AP Position (% Embryo Length)')
ylabel('Transition')
title('Transition from State0 to State1')

saveas(rfit_fig,['E:\Annika\mse_dissection\fig\mse_comparison_lateralML\exploratory_analyses\Hmm\rfit2_mean_fig.png'])

