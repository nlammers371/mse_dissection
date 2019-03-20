%generating a heatmap of fluorescence over space and time for each set ID
%Run from inside the exploratory_analyses folder
clear
close all
project = 'mse_comparison_mcp2f';
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
setID_vec = [nucleus_struct.setID];
setID_index = unique(setID_vec,'stable');
maxTime = ceil(max([nucleus_struct.time])/60); %finds the longest time in the 
%note: ceil rounds up to next integer - works well for binning
% make a reference vector indicating time and ap bins we want averages for
time_ref_vec = 1:maxTime;
ap_ref_vec = 1:100;
for s = 1:numel(setID_index) %repeat for each set
    sID = num2str(setID_index(s));
    setFilter = setID_vec==setID_index(s);
    fluo_vec = [nucleus_struct(setFilter).fluo];
    fluo_vec_zeros = fluo_vec;
    fluo_vec_zeros(isnan(fluo_vec_zeros)) = 0;
    time_vec = [nucleus_struct(setFilter).time]/60;
    ap_vec = [nucleus_struct(setFilter).ap_vector]*100;
    gName = unique({nucleus_struct(setFilter).genotype}, 'stable');
    if numel(gName)>1
        error 'Filtering error'
    end
    gNameStr = gName{1};
    %create two arrays which will hold our spot and nuclei intensity values
    spot_heatmap_array = NaN(numel(time_ref_vec),numel(ap_ref_vec));
    nuc_heatmap_array = NaN(numel(time_ref_vec),numel(ap_ref_vec));
    for a = 1:numel(ap_ref_vec)
        AP = ap_ref_vec(a);
        apFilter = round(ap_vec)==AP; % same idea: any entry that matches will have a 1
        for t = 1:maxTime
            time = time_ref_vec(t);
            tFilter = round(time_vec)==time;
            % now filter for relevant spot intensities using filter vectors
            fluo = fluo_vec(apFilter&tFilter);
            fluo_zeros = fluo_vec_zeros(apFilter&tFilter);
            spot_heatmap_array(t,a) = nanmean(fluo);
            nuc_heatmap_array(t,a) = mean(fluo_zeros);
        end
    end
    %now to make the heatmaps
    %first the average spot intensity
    spot_heatmap_fig = figure;
    imagesc(spot_heatmap_array)
    colorbar
    caxis([0,20e+5]); %make all the colorbars the same
    xlabel('AP position (% embryo length)')
    ylabel('Time (min)')
    c1 = colorbar;
    c1.Label.String = 'Average Spot Intensity (au)';
    title([gNameStr ' Spot Intensity Profile Set #' sID])
    saveas (spot_heatmap_fig,[figPath, gNameStr, '_nuc_heatmap_set', sID, '.png'])
    %then the average intensity across all nuclei
    nuc_heatmap_fig = figure;
    imagesc(nuc_heatmap_array)
    colorbar
    caxis([0,6e+5])
    %caxis to make all the colorbars the same
    xlabel('AP position (% embryo length)')
    ylabel('Time (min)')
    c2 = colorbar;
    c2.Label.String = 'Average Nuclear Intensity (au)';
    title([gNameStr, ' Nuclear Intensity Profile Set #', sID])
    saveas (nuc_heatmap_fig,[figPath, gNameStr, '_nuc_heatmap_set', sID, '.png']);
end
