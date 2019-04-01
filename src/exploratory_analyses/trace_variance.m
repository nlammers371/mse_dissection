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
for n = 1:numel(qc_nucleus_struct)
    timeVec(n,:) = qc_nucleus_struct(n).time_interp; %structure format (nucleus:time steps)
    apVec(n) = (qc_nucleus_struct(n).apMean)*100; %used mean AP so traces are only counted in one AP bin
    fluoVec(n,:) = qc_nucleus_struct(n).fluo_interp;
    gtypeVec(n) = {qc_nucleus_struct(n).genotype};
end

tFilterEarly = timeVec(1,:) < 1200; %"early" time points happen before 20 min
tFilterMid = timeVec(1,:) >= 1200 & timeVec(1,:) < 2400; %"mid" time points between 20 min and 40 min
tFilterLate = timeVec(1,:) >= 2400; %"late" time points are after 40 min
gIndex = {'Wt', 'Hb', 'Gt',};

for g = 1:numel(gIndex) %iterate over genotypes
    gFilter = strcmpi(gtypeVec, gIndex{g});
    fluo = fluoVec(gFilter,:);
    AP = apVec(gFilter);
    varFluo = nan(1,size(fluo,1));
    %make general trace variance plot over AP
    for i = 1:size(fluo,1)
        varFluo(i) = nanstd(fluo(i,:))/nanmean(fluo(i,:));
    end
    trace_variance_fig = figure;
    scatter(AP, varFluo)
    xlabel('AP position (%embryo length)')
    ylabel('Fluorescence Variance (standard deviation/mean)')
    title ([gIndex{g} ' Trace Variance of Active Nuclei'])
    grid on
    ylim([0,10])
    xlim([20,60])
    saveas (trace_variance_fig,[figPath gIndex{g} '_trace_variance_fig.png']);
    
    %Divide into time sections
    %Early
    eFluo = fluoVec(gFilter,tFilterEarly);
    eVarFluo = nan(1,size(eFluo,1));
    for i = 1:size(eFluo,1)
        eVarFluo(i) = nanstd(eFluo(i,:))/nanmean(eFluo(i,:));
    end
    etrace_variance_fig = figure;
    scatter(AP, eVarFluo)%,'FaceAlpha',.2)
    xlabel('AP position (%embryo length)')
    ylabel('Fluorescence Variance (standard deviation/mean)')
    title ([gIndex{g} ' Early Trace Variance of Active Nuclei (0 to 20 min)'])
    grid on
    ylim([0,10])
    xlim([20,60])
    saveas (etrace_variance_fig,[figPath gIndex{g} '_etrace_variance_fig.png']);
    
    %Middle
    mFluo = fluoVec(gFilter,tFilterMid);
    mVarFluo = nan(1,size(mFluo,1));
    for i = 1:size(mFluo,1)
        mVarFluo(i) = nanstd(mFluo(i,:))/nanmean(mFluo(i,:));
    end
    mtrace_variance_fig = figure;
    scatter(AP, mVarFluo)
    xlabel('AP position (%embryo length)')
    ylabel('Fluorescence Variance (standard deviation/mean)')
    title ([gIndex{g} ' Middle Trace Variance of Active Nuclei (20 to 40 min)'])
    grid on
    ylim([0,10])
    xlim([20,60])
    saveas (mtrace_variance_fig,[figPath gIndex{g} '_mtrace_variance_fig.png']);
    
    %Late
    lFluo = fluoVec(gFilter,tFilterLate);
    lVarFluo = nan(1,size(lFluo,1));
    for i = 1:size(lFluo,1)
        lVarFluo(i) = nanstd(lFluo(i,:))/nanmean(lFluo(i,:));
    end
    ltrace_variance_fig = figure;
    scatter(AP, lVarFluo)
    xlabel('AP position (%embryo length)')
    ylabel('Fluorescence Variance (standard deviation/mean)')
    title ([gIndex{g} ' Late Trace Variance of Active Nuclei (40 to 60 min)'])
    grid on
    ylim([0,10])
    xlim([20,60])
    saveas (ltrace_variance_fig,[figPath gIndex{g} '_ltrace_variance_fig.png']);
    
end

