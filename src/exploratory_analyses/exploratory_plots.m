%run from inside the exploratory_analyses folder
clear
close all
% specify project
project = 'mse_comparison_lateralML';
% this finds name of current directory
currentFolder = pwd;
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
% make a folder of same name in "fig" directory to hold any figures we make
figPath = ['../../fig/' project '/' fName '/']; % "../" tells code to go up one level
mkdir(figPath);
% specify path to data
dataPath = ['../../dat/' project '/'];
% load data
load([dataPath 'nucleus_struct.mat'])

% Now make some vectors with useful variables. The command
% [nucleus_struct.variable] will concatenate all variable values from
% across nuclei into one, long vector. This format is often easier to work
% with
fluo_vec = [nucleus_struct.fluo]; % spot intensity
fluo_vec_zeros = fluo_vec;
fluo_vec_zeros(isnan(fluo_vec_zeros)) = 0;
ap_vec = [nucleus_struct.ap_vector]*100; % AP position 
time_vec = [nucleus_struct.time]; % time
gtype_vec = [nucleus_struct.gtypeID]; % indicates genotype
% make index vectors to keep track of each unique genotype in the data
genotype_id_index = unique(gtype_vec,'stable')
genotype_str_index = unique({nucleus_struct.genotype},'stable') % if I don't include a semicolon it prints to command line

% First, lets just get a sense for the overall distribution of each
% variable of interest
% make and save a figure
fluo_hist = figure;
% plot
histogram(fluo_vec);
% add labels
xlabel('fluorescence (au)')
ylabel('count')
title('Distribution of fluorescence values')
% save
saveas(fluo_hist,[figPath 'fluo_hist.png'])

% same for AP position
ap_hist = figure;
% plot
histogram(ap_vec);
% add labels
xlabel('AP position (% embryo length)')
ylabel('count')
title('Distribution of AP positions')
% save
saveas(fluo_hist,[figPath 'ap_hist.png'])

% Typically, we'll probably want to break things up by genotype. We can use
% the index vectors to do this. Let's start by looking at average
% fluorescence by AP position

% make a reference vector indicating AP bins we want averages for
ap_ref_vec = 1:100;
% initialize an array to store fluo values for each AP-gtype combo
gtype_ap_fluo_array = NaN(numel(ap_ref_vec),numel(genotype_id_index));
% now we need to write for loops to calculate averages

% outer level: iterate over genotype
for g = 1:numel(genotype_id_index)
    gID = genotype_id_index(g);
    gName = genotype_str_index{g}; % note we use brackets to idex a string array
    gFilter = gtype_vec==gID; % vector of 1's and 0's. 1's indicate entry with correct genotype
    % inner level: iterate over AP
    for a = 1:numel(ap_ref_vec)
        AP = ap_ref_vec(a);
        apFilter = round(ap_vec)==AP; % same idea: any entry that matches will have a 1
        % now filter for relevant spot intensities using filters
        fluo = fluo_vec(gFilter&apFilter);
        fluo_zeros = fluo_vec_zeros(gFilter&apFilter);
        % record in array
        gtype_ap_fluo_array(a,g) = nanmean(fluo);
        gtype_ap_fluo_zeros_array(a,g) = mean(fluo_zeros);
    end
end

% make a figure for the average spot intensity (no zeros)
mean_fluo_fig = figure;
plot(ap_ref_vec,gtype_ap_fluo_array)
% make labels
xlabel('AP position (% embryo length)')
ylabel('spot intensity (au)')
title('Spatial spot intensity profiles by genotype')
% make a legend
legend(genotype_str_index{:})
% add grid lines
grid on
% save
saveas(mean_fluo_fig,[figPath 'mean_spot_fluo_fig.png'])

% make a figure for average nuclear intensity (with zeros)
mean_zeros_fluo_fig = figure;
plot(ap_ref_vec,gtype_ap_fluo_zeros_array)
% make labels
xlabel('AP position (% embryo length)')
ylabel('Nulcear intensity (au)')
title('Spatial Nuclear intensity profiles by genotype')
% make a legend
legend(genotype_str_index{:})
% add grid lines
grid on
% save
saveas(mean_zeros_fluo_fig,[figPath 'mean_nuc_fluo_fig.png'])