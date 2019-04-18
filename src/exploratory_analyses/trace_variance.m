%%need to fix this!!!
clear
close all
project = 'mse_comparison_lateralML';
currentFolder = pwd;
%find file name and make another file of the same name in the data folder
%to hold our figures
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
figPath = ['../../fig/' project '/' fName '/Variance/'];
mkdir(figPath);
%specify where the data is coming from and load the data
dataPath = ['../../dat/' project '/'];
load([dataPath 'qc_nucleus_struct.mat']);

%apply quality control filtering
qcFilter = [nucleus_struct.qc] == 1;
qc_nucleus_struct = nucleus_struct(qcFilter);


%create reference arrays
ap_ref_vec = 1:100;
gID_index = unique([qc_nucleus_struct.gtypeID],'stable');
gType_index = unique({qc_nucleus_struct.genotype}, 'stable');

varFluo = nan(1,numel(qc_nucleus_struct));
eVarFluo = nan(1,numel(qc_nucleus_struct));
mVarFluo = nan(1,numel(qc_nucleus_struct));
lVarFluo = nan(1,numel(qc_nucleus_struct));
apVec = nan(1,numel(qc_nucleus_struct));
gtypeVec = nan(1,numel(qc_nucleus_struct));

for n = 1:numel(qc_nucleus_struct)
    timeVec = [qc_nucleus_struct(n).time_interp]; 
    fluoVec = [qc_nucleus_struct(n).fluo_interp];
    gtypeVec(n) = mean(qc_nucleus_struct(n).gtypeID);
    apVec(n) = qc_nucleus_struct(n).apMean*100;
    
    varFluo(n) = nanstd(fluoVec)/nanmean(fluoVec);
    
    tFilterEarly = timeVec < 1200; %"early" time points happen before 20 min
    eFluo = fluoVec(tFilterEarly);
    eVarFluo(n) = nanstd(eFluo)/nanmean(eFluo);
    tFilterMid = timeVec >= 1200 & timeVec(1,:) < 2400; %"mid" time points between 20 min and 40 min
    mFluo = fluoVec(tFilterMid);
    mVarFluo(n) = nanstd(mFluo)/nanmean(mFluo);
    tFilterLate = timeVec >= 2400; %"late" time points are after 40 min
    lFluo = fluoVec(tFilterLate);
    lVarFluo(n) = nanstd(lFluo)/nanmean(lFluo);
end

meanVarFluo = nan(numel(gID_index), numel(ap_ref_vec));
SEMVarFluo = nan(numel(gID_index), numel(ap_ref_vec));

meaneVarFluo = nan(numel(gID_index),numel(ap_ref_vec));
SEMeVarFluo = nan(numel(gID_index),numel(ap_ref_vec));

meanmVarFluo = nan(numel(gID_index),numel(ap_ref_vec));
SEMmVarFluo = nan(numel(gID_index),numel(ap_ref_vec));

meanlVarFluo = nan(numel(gID_index),numel(ap_ref_vec));
SEMlVarFluo = nan(numel(gID_index),numel(ap_ref_vec));

for g = 1:numel(gID_index) %iterate over genotypes
    gFilter = gtypeVec == gID_index(g);
    gAPVec = apVec(gFilter);
    gVarFluo = varFluo(gFilter);
    geVarFluo = eVarFluo(gFilter);
    gmVarFluo = mVarFluo(gFilter);
    glVarFluo = lVarFluo(gFilter);
    
    %find averages over AP
    for j = 1:numel(ap_ref_vec)
        step = ap_ref_vec(j);
        stepFilter = round(gAPVec) == step;
        meanVarFluo(g,j) = nanmean(gVarFluo(stepFilter));
        SEMVarFluo(g,j) = nanstd(gVarFluo(stepFilter))/sqrt(sum(~isnan(gVarFluo(stepFilter))));
        meaneVarFluo(g,j) = nanmean(geVarFluo(stepFilter));
        SEMeVarFluo(g,j) = nanstd(geVarFluo(stepFilter))/sqrt(sum(~isnan(geVarFluo(stepFilter))));
        meanmVarFluo(g,j) = nanmean(gmVarFluo(stepFilter));
        SEMmVarFluo(g,j) = nanstd(gmVarFluo(stepFilter))/sqrt(sum(~isnan(gmVarFluo(stepFilter))));
        meanlVarFluo(g,j) = nanmean(glVarFluo(stepFilter));
        SEMlVarFluo(g,j) = nanstd(glVarFluo(stepFilter))/sqrt(sum(~isnan(glVarFluo(stepFilter))));
    end
    
    %plot broad trace variance patterns
    trace_variance_fig = figure;
    scatter(gAPVec, gVarFluo, 'MarkerEdgeAlpha', 0.5)
    xlabel('AP position (%embryo length)')
    ylabel('Trace Variance (standard deviation/mean)')
    title ([gType_index{g} ' Trace Variance of Active Nuclei'])
    grid on
    hold on
    %add the average variance over AP
    errorbar(ap_ref_vec, meanVarFluo(g,:),SEMVarFluo(g,:), 'LineWidth', 1.25)
    hold off
    xlim([20,50])
    ylim([0,1.3])
    saveas (trace_variance_fig,[figPath gType_index{g} '_trace_variance_fig.png']);
    
    %Do the same for the time intervals
    %Early
    etrace_variance_fig = figure;
    scatter(gAPVec, geVarFluo, 'MarkerEdgeAlpha', 0.5)
    xlabel('AP position (%embryo length)')
    ylabel('Fluorescence Variance (standard deviation/mean)')
    title ([gType_index{g} ' Early Trace Variance of Active Nuclei (0 to 20 min)'])
    grid on
    hold on
    %add the average variance over AP
    errorbar(ap_ref_vec, meaneVarFluo(g,:), SEMeVarFluo(g,:), 'LineWidth', 1.25)
    xlim([20,50])
    ylim([0,1.3])
    hold off
    saveas (etrace_variance_fig,[figPath gType_index{g} '_etrace_variance_fig.png']);
    
    %Middle
    mtrace_variance_fig = figure;
    scatter(gAPVec, gmVarFluo, 'MarkerEdgeAlpha', 0.5)
    xlabel('AP position (%embryo length)')
    ylabel('Fluorescence Variance (standard deviation/mean)')
    title ([gType_index{g} ' Middle Trace Variance of Active Nuclei (20 to 40 min)'])
    grid on
    hold on
    %add the average variance over AP
    errorbar(ap_ref_vec, meanmVarFluo(g,:), SEMmVarFluo(g,:), 'LineWidth', 1.25)
    xlim([20,50])
    ylim([0,1.3])
    hold off
    saveas (mtrace_variance_fig,[figPath gType_index{g} '_mtrace_variance_fig.png']);
    
    %Late
    ltrace_variance_fig = figure;
    scatter(gAPVec, glVarFluo, 'MarkerEdgeAlpha', 0.5)
    xlabel('AP position (%embryo length)')
    ylabel('Fluorescence Variance (standard deviation/mean)')
    title ([gType_index{g} ' Late Trace Variance of Active Nuclei (40 to 60 min)'])
    grid on
    hold on
    %add the average variance over AP
    errorbar(ap_ref_vec, meanlVarFluo(g,:), SEMlVarFluo(g,:), 'LineWidth', 1.25)
    xlim([20,50])
    ylim([0,1.3])
    hold off
    saveas (ltrace_variance_fig,[figPath gType_index{g} '_ltrace_variance_fig.png']);
end
%%
all_trace_var = figure;
for i = 1:numel(gID_index)
    errorbar(ap_ref_vec, meanVarFluo(i,:), SEMVarFluo(i,:), 'LineWidth', 1.25)
    hold on
end
hold off
xlabel('AP position (%embryo length)')
ylabel('Fluorescence Variance (standard deviation/mean)')
title ('Trace Variance of Active Nuclei')
grid on 
xlim([20,50])
ylim([0,1.3])
saveas (all_trace_var,[figPath 'all_trace_variance_fig.png']);
