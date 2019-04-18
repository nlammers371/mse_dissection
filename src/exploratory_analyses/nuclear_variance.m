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

%reference vectors
timeRef = 0:5:60; %specifies five minute time bins
apRef = 0:1:100; %specifies increments of 1%

%apply quality control filtering
qcFilter = [nucleus_struct.qc] == 1;
qc_nucleus_struct = nucleus_struct(qcFilter);

%create reference arrays
gID_index = unique([qc_nucleus_struct.gtypeID],'stable');
gType_index = unique({qc_nucleus_struct.genotype}, 'stable');
setID_index = unique([qc_nucleus_struct.setID], 'stable');
timeVec = [];
apVec = [];
fluoVec = [];
gtypeVec = [];
setIDVec = [];
for n = 1:numel(qc_nucleus_struct)
    timeVec = [timeVec (qc_nucleus_struct(n).time_interp)/60]; %want in terms of minutes
    apVec = [apVec (qc_nucleus_struct(n).ap_vector_interp)*100];
    fluoVec = [fluoVec qc_nucleus_struct(n).fluo_interp];
    gtypeVec(numel(gtypeVec)+1:numel(gtypeVec)+numel(qc_nucleus_struct(n).time_interp)) = mean(qc_nucleus_struct(n).gtypeID);
    setIDVec(numel(setIDVec)+1:numel(setIDVec)+numel(qc_nucleus_struct(n).time_interp)) = qc_nucleus_struct(n).setID;
end

%%
apVariance = NaN(numel(gType_index),numel(apRef));
apVariance0 = NaN(numel(gType_index),numel(apRef));

for g = 1:numel(gType_index) %iterate over genotypes
    variance_array = NaN(numel(timeRef),numel(apRef));
    variance_array0 = NaN(numel(timeRef),numel(apRef));
    gFilter = gtypeVec == gID_index(g);
    for a = 1:numel(apRef) %iterate over ap bins
        AP = apRef(a);
        apFilter = round(apVec)==AP;
        for t = 1:numel(timeRef) %iterate over time bins
            time = timeRef(t);
            tFilter = timeVec <= time & timeVec > (time-5); %300 sec = 5 min increments
            fluo = fluoVec(gFilter&apFilter&tFilter);
            variance_array(t,a) = nanstd(fluo,0,'all')/nanmean(fluo);
            fluo(isnan(fluo))=0;
            variance_array0(t,a) = nanstd(fluo,0,'all')/nanmean(fluo);
        end
        apFluo = fluoVec(gFilter&apFilter);
        apVariance(g,a) = (nanstd(apFluo,0,'all')/nanmean(apFluo));
        apFluo(isnan(apFluo))=0;
        apVariance0(g,a) = (std(apFluo,0,'all')/mean(apFluo));
    end
    variance_array_fig = figure;
    imagesc(variance_array)
    colorbar
    caxis([0,3])
    xlabel('AP position (% embryo length)')
    ylabel('Time (min)')
    title([gType_index{g}, ' Variance Across Active Nuclei'])
    saveas (variance_array_fig, [figPath, gType_index{g}, '_spot_variance_heatmap.png']);
    
    variance_array0_fig = figure;
    imagesc(variance_array0)
    colorbar
    caxis([0,7.5])
    xlabel('AP position (% embryo length)')
    ylabel('Time (min)')
    title([gType_index{g}, ' Variance Across All Nuclei'])
    saveas (variance_array0_fig, [figPath, gType_index{g}, '_nuc_variance_heatmap.png']);
end

ap_variance_fig = figure;
plot(apRef, apVariance, 'LineWidth', 1.25)
xlabel('AP position (%embryo length)')
ylabel('Variance (standard deviation/mean)')
title ('Fluorescence Variance Across Active Nuclei')
legend(gType_index{:})
grid on
saveas (ap_variance_fig,[figPath 'ap_spot_variance_fig.png']);

ap_variance0_fig = figure;
plot(apRef, apVariance0, 'LineWidth', 1.25)
xlabel('AP position (%embryo length)')
ylabel('Variance (standard deviation/mean)')
title ('Fluorescence Variance Across All Nuclei')
legend(gType_index{:})
grid on
saveas (ap_variance0_fig,[figPath 'ap_nuc_variance_fig.png']);

%%
%looking at variance with SEM calculated between sets of the same genotype
apVariance = NaN(numel(gType_index),numel(apRef));
apVarianceSEM = NaN(numel(gType_index),numel(apRef));

for g = 1:numel(gID_index)
    gFilter = gtypeVec == gID_index(g);
    for a = 1:numel(apRef)
        apFilter = round(apVec) == apRef(a);
        varianceVec = nan(1, numel(setID_index));
        for s = 1:numel(setID_index)
            setFilter = setIDVec == setID_index(s);
            fluo = fluoVec(gFilter&apFilter&setFilter);
            varianceVec(s) = nanstd(fluo)/nanmean(fluo);
        end
        apVariance(g,a) = nanmean(varianceVec);
        apVarianceSEM(g,a) = nanstd(varianceVec)/sqrt(sum(~isnan(varianceVec)));
    end
end

set_variance_fig = figure;
for i = 1:numel(gID_index)
   errorbar(apRef, apVariance(i,:),apVarianceSEM(i,:), 'LineWidth', 1.25) 
   hold on
end
hold off
xlabel('AP position (%embryo length)')
ylabel('Variance (standard deviation/mean)')
title ('Fluorescence Variance Between Active Nuclei')
legend(gType_index{:})
grid on
saveas (set_variance_fig,[figPath 'set_nuc_variance_fig.png']);
    