%Stripe Center defined for this analysis as being between .38 and .42
%(fraction of embryo length)

%run from inside exploratory analysis folder
clear
close all
project = 'mse_comparison_lateralML';
currentFolder = pwd;
%find file name and make another file of the same name in the data folder
%to hold our figures
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
figPath = ['../../fig/' project '/' fName '/'];
mkdir(figPath);
%specify where the data is coming from and load the data
dataPath = ['../../dat/' project '/'];
load([dataPath 'qc_nucleus_struct.mat']);

ncFilter = [nucleus_struct.nc_qc_flag] == 1;
qc_nucleus_struct = nucleus_struct(ncFilter);

%%
%define useful vectors
fluo_vec = [qc_nucleus_struct.fluo_interp];
fluo_vec_zeros = fluo_vec;
fluo_vec_zeros(isnan(fluo_vec_zeros)) = 0;
time_vec = [qc_nucleus_struct.time_interp]/60;
ap_vec = [qc_nucleus_struct.ap_vector_interp]*100;
gtype_vec = [qc_nucleus_struct.gtypeID_interp];
%index vectors to keep track of individual genotypes
gID_index = unique(gtype_vec,'stable');
genotype_str_index = unique({qc_nucleus_struct.genotype}, 'stable');

apFilter = ap_vec > 38 & ap_vec < 42;
maxTime = 60; 
time_ref_vec = 1:maxTime;

%make an array to hold all of our average fluorescences over time per
%genotype
avespotfluo_array = NaN(numel(time_ref_vec),numel(gID_index)); %array for average intensity of "on" spots
avenucfluo_array = NaN(numel(time_ref_vec),numel(gID_index)); %array for average intensity for all nuclei
totfluo_array = NaN(numel(time_ref_vec),numel(gID_index)); %array for totals

avespotfluo_SEM = NaN(numel(time_ref_vec),numel(gID_index)); 
avenucfluo_SEM = NaN(numel(time_ref_vec),numel(gID_index)); 

%iteratively find the averages per genotype for each bin
for g = 1:numel(gID_index)%iterate over genotypes
    gID = gID_index(g);
    gName = genotype_str_index{g};
    gFilter = gtype_vec==gID; 
    for t = 1:numel(time_ref_vec) %iterate over time
        time = time_ref_vec(t);
        tFilter = round(time_vec)==time;
        fluo = fluo_vec(gFilter&tFilter&apFilter); %filter for relevant spot intensities
        fluo_zeros = fluo_vec_zeros(gFilter&tFilter&apFilter);
        %record the averages in the array
        avespotfluo_array(t,g) = nanmean(fluo); %average without NaNs
        avespotfluo_SEM(t,g) = nanstd(fluo)/sqrt(sum(~isnan(fluo)));
        avenucfluo_array(t,g) = mean(fluo_zeros); %average with zeros
        avenucfluo_SEM(t,g) = std(fluo_zeros)/sqrt(numel(fluo_zeros));
        totfluo_array(t,g) = nansum(fluo); %add without NaNs
    end
end
%plot the average spot intensity as a function of time
meanspotfluo_time_fig = figure;%make figure
for i=1:numel(gID_index)
    errorbar(time_ref_vec, avespotfluo_array(:,i), avespotfluo_SEM(:,i), 'LineWidth', 1.25)
    hold on
end
hold off
xlabel('Time (minutes)')
ylabel('Average Spot Intensity (au)')
title ('Stripe Proper Spot Intensity by Genotype')
legend(genotype_str_index{:})
grid on
saveas (meanspotfluo_time_fig,[figPath 'SPmean_spot_fluo_time_fig.png']);

%plot the total fluorescence as a function of time
% tot_fluo_time_fig = figure;
% for i=1:numel(gID_index)
%     errorbar(time_ref_vec, totfluo_array(:,i), totfluo_SEM(:,i))
%     hold on
% end
% hold off
% xlabel('Time (minutes)')
% ylabel('Total Spot Intensity (au)')
% title ('Stripe Center Total Intensity by Genotype')
% legend(genotype_str_index{:})
% grid on
%saveas (tot_fluo_time_fig,[figPath 'SCtot_fluo_time_fig.png']);

%plot the average nuclear fluorescence as a function of time
meannucfluo_time_fig = figure;%make figure
for i=1:numel(gID_index)
    errorbar(time_ref_vec, avenucfluo_array(:,i), avenucfluo_SEM(:,i), 'LineWidth', 1.25)
    hold on
end
hold off
xlabel('Time (minutes)')
ylabel('Average Nuclear Intensity (au)')
title ('Stripe Proper Nuclear Intensity by Genotype')
legend(genotype_str_index{:})
ylim([0,6e5])
grid on
saveas (meannucfluo_time_fig,[figPath 'SPmean_nuc_fluo_time_fig.png']);

%%
%%%look at the fraction of active nuclei over time
ncFilter = [nucleus_struct.nc_qc_flag] == 1;
nc_nucleus_struct = nucleus_struct(ncFilter);

onTimeVec = nan(1,numel(nc_nucleus_struct));
offTimeVec = nan(1,numel(nc_nucleus_struct));
meanAPVec = nan(1,numel(nc_nucleus_struct));
gIDVec = nan(1,numel(nc_nucleus_struct));

gID_index = unique([nc_nucleus_struct.gtypeID],'stable');
genotype_str_index = unique({nucleus_struct.genotype}, 'stable');
timeRefVec = 1:60;


for n = 1:numel(nc_nucleus_struct)
    if sum(~isnan(nc_nucleus_struct(n).fluo_interp))>0
        onTimeVec(n) = nc_nucleus_struct(n).time_interp(find(~isnan(nc_nucleus_struct(n).fluo_interp),1))/60;
        offTimeVec(n) =  nc_nucleus_struct(n).time_interp(find(~isnan(nc_nucleus_struct(n).fluo_interp),1,'last'))/60;
    end
    meanAPVec(n) = nc_nucleus_struct(n).apMean*100;
    gIDVec(n) = mean(nc_nucleus_struct(n).gtypeID);
end

fracton_array = nan(numel(timeRefVec),numel(gID_index));
apFilter = meanAPVec > 38 & meanAPVec < 42;
for g=1:numel(gID_index)
    gFilter = gIDVec == gID_index(g);
    for t=1:numel(timeRefVec)
        tFilter = onTimeVec < timeRefVec(t) & offTimeVec > timeRefVec(t);
        fracton_array(t,g) = sum(tFilter&gFilter&apFilter)/sum(gFilter&apFilter);
    end
end

fract_on_fig = figure;
plot(timeRefVec,fracton_array)
hold off
xlabel('Time (minutes)')
ylabel('Fraction of Active Nuclei')
title ('Fraction of Active Nuclei over Time')
legend(genotype_str_index{:})
grid on
saveas (fract_on_fig,[figPath 'SPfract_on_fig.png']);

%%
multiple_array = nan(size(fracton_array,1),size(fracton_array,2));
for x = 1:size(fracton_array,1)
    for y = 1:size(fracton_array,2)
        multiple_array(x,y) = avespotfluo_array(x,y) * fracton_array(x,y);
    end
end
fig = figure;
plot(timeRefVec, multiple_array)
xlabel('Time (minutes)')
ylabel('Average Nuclear Intensity (au)')
title ('spot x active nuclei by Genotype')
legend(genotype_str_index{:})
grid on
