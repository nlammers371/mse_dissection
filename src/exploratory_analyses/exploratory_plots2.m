%Run from inside the exploratory_analyses folder
clear
close all
project = 'mse_comparison_mcp4f';
currentFolder = pwd;
%find file name and make another file of the same name in the data folder
%to hold our figures
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
figPath = ['../../fig/' project '/' fName '/'];
mkdir(figPath);
%specify where the data is coming from and load the data
dataPath = ['../../dat/' project '/'];
load([dataPath 'nucleus_struct.mat']);
%define useful vectors
fluo_vec = [nucleus_struct.fluo];
fluo_vec_zeros = fluo_vec;
fluo_vec_zeros(isnan(fluo_vec_zeros)) = 0;
time_vec = [nucleus_struct.time];
ap_vec = [nucleus_struct.ap_vector]*100;
gtype_vec = [nucleus_struct.gtypeID];
%index vectors to keep track of individual genotypes
genotype_id_index = unique(gtype_vec,'stable');
genotype_str_index = unique({nucleus_struct.genotype}, 'stable');
%stable keeps things in the same order as they originally appeared so that
%the index numbers match up with the string identifiers

%I want to know what the differences in spot timing are between the three
%genotypes (average fluoresence for each genotype over time)
%first need to find average and total fluoresence at each time "bin" we create
time_vec_minutes = time_vec/60; %time in terms of minutes
maxTime = ceil(max(time_vec_minutes)); %finds the longest time in the 
%note: ceil rounds up to next integer - works well for binning
% make a reference vector indicating time bins we want averages for
time_ref_vec = 1:maxTime;
%%
%make an array to hold all of our average fluorescences over time per
%genotype
gtype_time_avespotfluo_array = NaN(numel(time_ref_vec),numel(genotype_id_index)); %array for average intensity of "on" spots
gtype_time_avenucfluo_array = NaN(numel(time_ref_vec),numel(genotype_id_index)); %array for average intensity for all nuclei
gtype_time_totfluo_array = NaN(numel(time_ref_vec),numel(genotype_id_index)); %array for totals
%iteratively find the averages per genotype for each bin
for g = 1:numel(genotype_id_index)%iterate over genotypes
    gID = genotype_id_index(g);
    gName = genotype_str_index{g};
    gFilter = gtype_vec==gID; % vector of 1's and 0's. 1's indicate entry with correct genotype
    for t = 1:numel(time_ref_vec); %iterate over time
        time = time_ref_vec(t);
        tFilter = round(time_vec_minutes)==time;
        fluo = fluo_vec(gFilter&tFilter); %filter for relevant spot intensities
        fluo_zeros = fluo_vec_zeros(gFilter&tFilter);
        %record the averages in the array
        gtype_time_avespotfluo_array(t,g) = nanmean(fluo); %average without NaNs
        gtype_time_avenucfluo_array(t,g) = mean(fluo_zeros); %average with zeros
        gtype_time_totfluo_array(t,g) = nansum(fluo); %add without NaNs
    end
end
%plot the average spot intensity as a function of time
meanspotfluo_time_fig = figure;%make figure
plot(time_ref_vec,gtype_time_avespotfluo_array)
xlabel('Time (minutes)')
ylabel('Average Spot Intensity (au)')
title ('Temporal Spot Intensity Profiles by Genotype')
legend(genotype_str_index{:})
grid on
saveas (meanspotfluo_time_fig,[figPath 'mean_spot_fluo_time_fig.png']);
%plot the total fluorescence as a function of time
tot_fluo_time_fig = figure;
plot(time_ref_vec,gtype_time_totfluo_array)
xlabel('Time (minutes)')
ylabel('Total Spot Intensity (au)')
title ('Total Temporal Intensity Profiles by Genotype')
legend(genotype_str_index{:})
grid on
saveas (tot_fluo_time_fig,[figPath 'tot_fluo_time_fig.png']);
%plot the average nuclear fluorescence as a function of time
meannucfluo_time_fig = figure;%make figure
plot(time_ref_vec,gtype_time_avenucfluo_array)
xlabel('Time (minutes)')
ylabel('Average Nuclear Intensity (au)')
title ('Temporal Nuclear Intensity Profiles by Genotype')
legend(genotype_str_index{:})
grid on
saveas (meannucfluo_time_fig,[figPath 'mean_nuc_fluo_time_fig.png']);
%%
%since I just want to compare gt, wt, and hb spot fluorescence over time
gtype_time_avespotfluo_array = NaN(numel(time_ref_vec),3); %array for average intensity of "on" spots
gtype_time_avenucfluo_array = NaN(numel(time_ref_vec),3); %array for average intensity for all nuclei
for g = 1:3%iterate over these three genotypes
    gID = genotype_id_index(g);
    gName = genotype_str_index{g};
    gFilter = gtype_vec==gID; % vector of 1's and 0's. 1's indicate entry with correct genotype
    for t = 1:numel(time_ref_vec); %iterate over time
        time = time_ref_vec(t);
        tFilter = round(time_vec_minutes)==time;
        fluo = fluo_vec(gFilter&tFilter); %filter for relevant spot intensities
        fluo_zeros = fluo_vec_zeros(gFilter&tFilter);
        %record the averages in the array
        gtype_time_avespotfluo_array(t,g) = nanmean(fluo); %average without NaNs
        gtype_time_avenucfluo_array(t,g) = mean(fluo_zeros); %average with zeros
        gtype_time_totfluo_array(t,g) = nansum(fluo); %add without NaNs
    end
end
%plot the average spot intensity as a function of time
meanspotfluo_time_fig = figure;%make figure
plot(time_ref_vec,gtype_time_avespotfluo_array)
xlabel('Time (minutes)')
ylabel('Average Spot Intensity (au)')
title ('Temporal Spot Intensity Profiles by Genotype')
legend(genotype_str_index{:})
grid on
saveas (meanspotfluo_time_fig,[figPath 'mean_spot_fluo_time_fig2.png']);
