%generating a scatter plot of nuclear on times across the Ap axis for each genotype
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

%First we need vectors of all the nuclear on times, ap positions, and genotypes
nuc_on_time_vec = NaN(1,numel(nucleus_struct));
nuc_genotype_vec = NaN(1,numel(nucleus_struct));
nuc_ap_vec = NaN(1,numel(nucleus_struct));
nuc_setID_vec = NaN(1,numel(nucleus_struct));
for n = 1:numel(nucleus_struct)
    fluo = nucleus_struct(n).fluo;
    time = nucleus_struct(n).time;
    first_on = find(~isnan(fluo),1);
    if ~isempty(first_on)
        on_time = time(first_on)/60; %convert to minutes
    else 
        on_time = NaN;
    end
    nuc_on_time_vec(n) = on_time;
    nuc_genotype_vec(n) = mean(nucleus_struct(n).gtypeID); %took average because there are multiple entries for the gtypeID
    nuc_ap_vec(n) = mean(nucleus_struct(n).ap_vector)*100; %average AP position as a percent of embryo length
    nuc_setID_vec(n) = nucleus_struct(n).setID;
end
%%
%now we can make some histograms for each genotype
for g = 1:numel(genotype_id_index)
    gID = genotype_id_index(g);
    gName = genotype_str_index{g};
    %need to  filter out everything to just focus on one genotype at a time
    gFilter = nuc_genotype_vec==gID; % vector of 1's and 0's. 1's indicate entry with correct genotype
    
    fig_nuc_on_time = nuc_on_time_vec(gFilter);
    on_time_hist = figure;
    histogram(fig_nuc_on_time)
    title([gName, 'Distribution of Nuclear On-Times']);
    xlabel('Time (min)')
    ylabel('Count')
    saveas(on_time_hist,[figPath gName '_on_time_hist.png'])
end

%%
%Or we can make scatter plots over AP and time where each point = 1 nucleus
for g = 1:numel(genotype_id_index)
    gID = genotype_id_index(g);
    gName = genotype_str_index{g};
    %need to  filter out everything to just focus on one genotype at a time
    gFilter = nuc_genotype_vec==gID;
    gtype_nuc_on_time = nuc_on_time_vec(gFilter);
    gtype_nuc_ap = nuc_ap_vec(gFilter);
    ap_ontime_plot = figure;
    scatter(gtype_nuc_ap,gtype_nuc_on_time, 'filled');
    xlim([15,65]);
    ylim([0,60]);
    title([gName, ' Nuclear On Times Across AP']);
    xlabel('AP position (% embryo length)');
    ylabel('Nuclear On-Time (min)');
    saveas(ap_ontime_plot,[figPath gName '_ap_ontime_plot.png'])
end
    



%%
%If I want to make scatter plots with different colors for each data set
for g = 1:numel(genotype_id_index)
    gID = genotype_id_index(g);
    gName = genotype_str_index{g};
    %need to  filter out everything to just focus on one genotype at a time
    gFilter = nuc_genotype_vec==gID;
    gtype_nuc_on_time = nuc_on_time_vec(gFilter);
    gtype_nuc_ap = nuc_ap_vec(gFilter);
    gtype_nuc_setID = nuc_setID_vec(gFilter);
    gtype_nuc_setID_index = unique(gtype_nuc_setID,'stable');
    if numel(gtype_nuc_setID_index) > 1
        %make a structure to allow for group scatter plots
        on_time_struct = struct('OnTime',{},'ApPosition',{},'SetID',{});
        for a = 1:numel(gtype_nuc_ap)
            %fill the structure
            on_time_struct(a).OnTime = gtype_nuc_on_time(a);
            on_time_struct(a).ApPosition = gtype_nuc_ap(a);
            on_time_struct(a).SetID = gtype_nuc_setID(a);
        end
        %make a plot separated by set number
        ap_ontime_set_plot = figure;
        gscatter([on_time_struct.ApPosition], [on_time_struct.OnTime],[on_time_struct.SetID], 'bgrk', '.', 10);
        xlim([15,65]);
        ylim([0,60]);
        title([gName, ' Nuclear On Times Across AP']);
        xlabel('AP position (% embryo length)');
        ylabel('Nuclear On-Time (min)');
        saveas(ap_ontime_set_plot,[figPath gName '_ap_ontime_plot_set.png'])
        %make a normal scatter plot of all the data combined
        ap_ontime_plot = figure;
        scatter(gtype_nuc_ap,gtype_nuc_on_time, 'filled');
        xlim([15,65]);
        ylim([0,60]);
        title([gName, ' Nuclear On Times Across AP']);
        xlabel('AP position (% embryo length)');
        ylabel('Nuclear On-Time (min)');
        saveas(ap_ontime_plot,[figPath gName '_ap_ontime_plot.png'])
    else
        ap_ontime_plot = figure;
        scatter(gtype_nuc_ap,gtype_nuc_on_time, 'filled');
        xlim([15,65]);
        ylim([0,60]);
        title([gName, ' Nuclear On Times Across AP']);
        xlabel('AP position (% embryo length)');
        ylabel('Nuclear On-Time (min)');
        saveas(ap_ontime_plot,[figPath gName '_ap_ontime_plot.png'])
    end
end


