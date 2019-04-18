%Run from inside the exploratory_analyses folder
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

%create reference indexes and vectors:
APRefVec = 1:100;
tRefVec = 1:60;

gID_index = unique([nucleus_struct.gtypeID],'stable');
gType_index = unique({nucleus_struct.genotype}, 'stable');
fractOn = NaN(numel(APRefVec), numel(gID_index));
gIDVec = NaN(1, numel(nucleus_struct));
meanAPVec = NaN(1, numel(nucleus_struct));
onVec = NaN(1, numel(nucleus_struct));

for i=1:numel(nucleus_struct)
    gIDVec(i) = mean(nucleus_struct(i).gtypeID);
    meanAPVec(i) = nucleus_struct(i).apMean*100;
    if sum(~isnan(nucleus_struct(i).fluo))>0 %if there is any recorded fluorescence
        onVec(i) = 1;
    end
end

for a=1:numel(APRefVec)
    apFilter = round(meanAPVec)==a;
    for g=1:numel(gType_index)
        gFilter = gIDVec == gID_index(g);
        nuclei = onVec(apFilter&gFilter);
        fractOn(a,g) = nansum(nuclei)/numel(nuclei);
    end
end

frac_on_fig = figure;
plot(APRefVec, fractOn, 'LineWidth',1.25)
xlabel('AP Position(%embryo length)')
ylabel('Fraction of Active Nuclei')
title('Fraction of Active Nuclei over AP')
legend(gType_index{:})
%saveas(frac_on_fig, [figPath 'faction_active.png']);

