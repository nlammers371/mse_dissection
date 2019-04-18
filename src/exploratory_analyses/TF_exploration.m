%Run from inside the exploratory_analyses folder
clear
close all
project = 'mse_comparison_lateralML';
currentFolder = pwd;
%find file name and make another file of the same name in the data folder
%to hold our figures
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
figPath = ['../../fig/' project '/' fName '/FluoBins/'];
mkdir(figPath);
%specify where the data is coming from and load the data
dataPath = ['../../dat/' project '/'];
load([dataPath 'nucleus_struct_protein.mat'])

nucFilter = [nucleus_struct.nc_qc_flag] == true;
qc_nucleus_struct = nucleus_struct(nucFilter);

apRefVec = 1:100;
meanGtVec = NaN(1,numel(apRefVec));
meanKrVec = NaN(1,numel(apRefVec));
meanHbVec = NaN(1,numel(apRefVec));
meanBcdVec = NaN(1,numel(apRefVec));
apVec = [];
GiantVec = [];
KruppelVec =  [];
HunchbackVec = [];
BicoidVec = [];

for n= 1:numel(qc_nucleus_struct)
    apVec = [apVec [qc_nucleus_struct(n).ap_vector_interp]*100];
    GiantVec = [GiantVec qc_nucleus_struct(n).tf_array(:,1).'];
    KruppelVec =  [KruppelVec qc_nucleus_struct(n).tf_array(:,2).'];
    HunchbackVec = [HunchbackVec qc_nucleus_struct(n).tf_array(:,3).'];
    BicoidVec = [BicoidVec qc_nucleus_struct(n).tf_array(:,4).'];
end

for a = 1:numel(apRefVec)
    apFilter = round(apVec) == apRefVec(a);
    Giant = GiantVec(apFilter);
    Kruppel = KruppelVec(apFilter);
    Hunchback = HunchbackVec(apFilter);
    Bicoid = BicoidVec(apFilter);
    
    meanGtVec(a) = nanmean(Giant);
    meanKrVec(a) = nanmean(Kruppel);
    meanHbVec(a) = nanmean(Hunchback);
    meanBcdVec(a) = nanmean(Bicoid);
end

tfFig = figure;
plot(apRefVec, meanGtVec)
hold on
plot(apRefVec, meanKrVec)
hold on
plot (apRefVec, meanHbVec)
hold on 
plot (apRefVec, meanBcdVec)
xlabel('Ap Position')