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
for n = 1:numel(qc_nucleus_struct)
    timeVec(n,:) = qc_nucleus_struct(n).time_interp; %structure format (nucleus:time steps)
    apVec(n,:) = (qc_nucleus_struct(n).ap_vector_interp)*100;
    fluoVec(n,:) = qc_nucleus_struct(n).fluo_interp;
    gtypeVec(n,1:181) = {qc_nucleus_struct(n).genotype};
end

gIndex = {'Wt', 'Hb', 'Gt'};
%%
apVariance = NaN(numel(gIndex),numel(apRef));
apVariance0 = NaN(numel(gIndex),numel(apRef));
for g = 1:numel(gIndex) %iterate over genotypes
    variance_array = NaN(numel(timeRef),numel(apRef));
    variance_array0 = NaN(numel(timeRef),numel(apRef));
    gFilter = strcmpi(gtypeVec, gIndex{g});
    for a = 1:numel(apRef) %iterate over ap bins
        AP = apRef(a);
        apFilter = round(apVec(:,:))==AP;
        for t = 1:numel(timeRef) %iterate over time bins
            time = timeRef(t)*60;
            tFilter = timeVec <= time & timeVec > (time-300); %300 sec = 5 min increments
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
    caxis([0,7.5])
    xlabel('AP position (% embryo length)')
    ylabel('Time (min)')
    title([gIndex{g}, ' Variance Across Active Nuclei'])
    saveas (variance_array_fig, [figPath, gIndex{g}, '_spot_variance_heatmap.png']);
    
    variance_array0_fig = figure;
    imagesc(variance_array0)
    colorbar
    caxis([0,7.5])
    xlabel('AP position (% embryo length)')
    ylabel('Time (min)')
    title([gIndex{g}, ' Variance Across All Nuclei'])
    saveas (variance_array0_fig, [figPath, gIndex{g}, '_nuc_variance_heatmap.png']);
end

ap_variance_fig = figure;
plot(apRef, apVariance)
xlabel('AP position (%embryo length)')
ylabel('Variance (standard deviation/mean)')
title ('Fluorescence Variance Across Active Nuclei')
legend(gIndex{:})
grid on
saveas (ap_variance_fig,[figPath 'ap_spot_variance_fig.png']);

ap_variance0_fig = figure;
plot(apRef, apVariance0)
xlabel('AP position (%embryo length)')
ylabel('Variance (standard deviation/mean)')
title ('Fluorescence Variance Across All Nuclei')
legend(gIndex{:})
grid on
saveas (ap_variance0_fig,[figPath 'ap_nuc_variance_fig.png']);