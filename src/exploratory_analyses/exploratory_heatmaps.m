%generating a heatmap of fluorescence over space and time
%Run from inside the exploratory_analyses folder
clear
close all
project = 'mse_comparison_mcp4f';
currentFolder = pwd;
%find file name and make another file of the same name in the data folder
%to hold our figures
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
figPath = ['../../fig/' project '/' fName '/Heatmaps/'];
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
time_vec_minutes = time_vec/60; %time in terms of minutes
maxTime = ceil(max(time_vec_minutes)); %finds the longest time in the 
%note: ceil rounds up to next integer - works well for binning
% make a reference vector indicating time bins we want averages for
time_ref_vec = 1:maxTime;
ap_ref_vec = 1:100;
%Need to iterate over genotype, AP position, and time
for g = 1:numel(genotype_id_index)
    %create two arrays which will hold our spot and nuclei intensity values
    spot_heatmap_array = NaN(numel(time_ref_vec),numel(ap_ref_vec));
    nuc_heatmap_array = NaN(numel(time_ref_vec),numel(ap_ref_vec));
    gID = genotype_id_index(g);
    gName = genotype_str_index{g};
    %need to  filter out everything to just focus on one genotype at a time
    gFilter = gtype_vec==gID; % vector of 1's and 0's. 1's indicate entry with correct genotype
    for a = 1:numel(ap_ref_vec)
        AP = ap_ref_vec(a);
        apFilter = round(ap_vec)==AP; % same idea: any entry that matches will have a 1
        for t = 1:maxTime
            time = time_ref_vec(t);
            tFilter = round(time_vec_minutes)==time;
            % now filter for relevant spot intensities using filter vectors
            fluo = fluo_vec(gFilter&apFilter&tFilter);
            fluo_zeros = fluo_vec_zeros(gFilter&apFilter&tFilter);
            spot_heatmap_array(t,a) = nanmean(fluo);
            nuc_heatmap_array(t,a) = mean(fluo_zeros);
        end
    end
    %make both figures for each genotype
    %first the average spot intensity
    spot_heatmap_fig = figure
    imagesc(spot_heatmap_array)
    colorbar
    caxis([0,20e+5]); %make all the colorbars the same
    xlabel('AP position (% embryo length)')
    ylabel('Time (min)')
    c1 = colorbar;
    c1.Label.String = 'Average Spot Intensity (au)';
    title([gName, ' Spatial Temporal Spot Intensity Profile'])
    saveas (spot_heatmap_fig,[figPath, gName, '_spot_heatmap_fig.png']);
    %then the average intensity across all nuclei
    nuc_heatmap_fig = figure
    imagesc(nuc_heatmap_array)
    colorbar
    caxis([0,6e+5])
    %caxis to make all the colorbars the same
    xlabel('AP position (% embryo length)')
    ylabel('Time (min)')
    c2 = colorbar;
    c2.Label.String = 'Average Nuclear Intensity (au)';
    title([gName, ' Spatial Temporal Nuclear Intensity Profile'])
    saveas (nuc_heatmap_fig,[figPath, gName, '_nuc_heatmap_fig.png']);
end
            
        
        

