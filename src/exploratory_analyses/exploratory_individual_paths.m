%generate a fluorescence over time trace for each individual nucleus that
%actually has a spot (discards nuclei with no spots)
%run from inside the exploratory analyses folder
clear 
close all
project = 'mse_comparison_lateralML';
currentFolder = pwd;
%specify folders for the data and for the figures we make
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
figPath = ['../../fig/' project '/' fName '/IndividualTraces/'];
mkdir(figPath);
dataPath = ['../../dat/' project '/'];
%load the data!
load([dataPath 'nucleus_struct.mat']);
%define and index genotype vectors
gtype_vec = [nucleus_struct.gtypeID];
time_vec = [nucleus_struct.time]/60; %divide by 60 to get time in min
fluo_vec = [nucleus_struct.fluo];
genotype_id_index = unique(gtype_vec,'stable');
genotype_str_index = unique({nucleus_struct.genotype}, 'stable');
%make genotype-specific folders
for g = 1:numel(genotype_id_index)
    gName = genotype_str_index{g};
    gtypePath = [figPath, '/' gName '/'];
    mkdir(gtypePath);
end
        
%%
%iterate over individual nuclei to create traces and sort to folders
for n=1:numel(nucleus_struct)
    fluo = [nucleus_struct(n).fluo];
    time = [nucleus_struct(n).time]/60;
    gID = mean(nucleus_struct(n).gtypeID); %took average because there are multiple entries for the gtypeID
    gName = genotype_str_index{gID};
    if sum(~isnan(fluo)) >= 10 %if the number of non-NaN values is more than 10, make a figure
        n_str = num2str(n);
        particle_trace = figure('visible','off');
        plot(time,fluo);
        title(['Particle Trace for ' gName ' Nucleus ' n_str]);
        xlabel('Time (min)');
        xlim([0,60]);
        ylabel('Intensity (au)');
        saveas(particle_trace,[figPath, '/' gName '/particle_trace_' n_str '.png']);
    end

end
