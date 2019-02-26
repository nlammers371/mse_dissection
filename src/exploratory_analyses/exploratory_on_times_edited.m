%generating a scatter plot of nuclear on times across the Ap axis for each
%genotype, color coded by set ID
%Run from inside the exploratory_analyses folder
clear
close all
project = 'mse_comparison_mcp4f';
currentFolder = pwd;
%find file name and make another file of the same name in the data folder
%to hold our figures
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
figPath = ['../../fig/' project '/' fName '/OnTimes/'];
mkdir(figPath);
%specify where the data is coming from and load the data
dataPath = ['../../dat/' project '/'];
load([dataPath 'nucleus_struct.mat']);
%define and index genotype vectors
gtype_vec = [nucleus_struct.gtypeID];
genotype_id_index = unique(gtype_vec,'stable');
genotype_str_index = unique({nucleus_struct.genotype}, 'stable');
%stable keeps things in the same order as they originally appeared so that
%the index numbers match up with the string identifiers

%now we need to make a structure to hold all the parameters we are going to
%extract for each individual nucleus
on_time_struct = struct('onTime',{},'apPosition',{},'setID',{},'gtypeID',{});
for n = 1:numel(nucleus_struct)
    fluo = nucleus_struct(n).fluo; %pulls out all the fluorescence values for the particular nucleus
    time = nucleus_struct(n).time; %same for time
    first_on = find(~isnan(fluo),1); %find the first entry that actually records a spot
    if ~isempty(first_on)
        on_time = time(first_on)/60;%convert to minutes
    else %if there are no spots recorded for the nucleus
        on_time = NaN;
    end
    on_time_struct(n).onTime = on_time;
    on_time_struct(n).apPosition = mean(nucleus_struct(n).ap_vector)*100; %average AP position as a percent of embryo length
    on_time_struct(n).setID = nucleus_struct(n).setID;
    on_time_struct(n).gtypeID = unique(nucleus_struct(n).gtypeID); %there should only be one unique gtypeID
end
%%
%now that we have a structure, we can start filtering by genotype and making
%scatter plots
for g = 1:numel(genotype_id_index)
    gID = genotype_id_index(g);
    gName = genotype_str_index{g};
    %need to  filter out everything to just focus on one genotype at a time
    gFilter = [on_time_struct.gtypeID] == gID;
    %make a scatter plot that is color-coded by set id 
    ap_ontime_set_plot = figure;
    gscatter([on_time_struct(gFilter).apPosition], [on_time_struct(gFilter).onTime],[on_time_struct(gFilter).setID], 'bgrk', '.', 10)
    xlim([15,65]);
    ylim([0,60]);
    title([gName, ' Nuclear On-Times Across AP']);
    xlabel('AP position (% embryo length)');
    ylabel('Nuclear On-Time (min)');
    %saveas(ap_ontime_set_plot,[figPath gName '_ap_ontime_plot_set.png'])
    %make a scatter plot of the overall trends (not color-coded_
    ap_ontime_plot = figure;
    scatter([on_time_struct(gFilter).apPosition],[on_time_struct(gFilter).onTime], 'filled');
    xlim([15,65]);
    ylim([0,60]);
    title([gName, ' Nuclear On Times Across AP']);
    xlabel('AP position (% embryo length)');
    ylabel('Nuclear On-Time (min)');
    %saveas(ap_ontime_plot,[figPath gName '_ap_ontime_plot.png'])
end    
