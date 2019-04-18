%run inside the exploratory analyses folder

clear
close all
project = 'mse_comparison_lateralML';
currentFolder = pwd;
%find file name and make another file of the same name in the data folder
%to hold our figures
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
figPath = ['../../fig/' project '/' fName '/QCHeatmaps/'];
mkdir(figPath);
%specify where the data is coming from and load the data
dataPath = ['../../dat/' project '/'];
load([dataPath 'qc_nucleus_struct.mat']);

%reference vectors
timeRef = 1:60;
apRef = 1:100;

%apply quality control filtering
qcFilter = [nucleus_struct.nc_qc_flag] == 1;
qc_nucleus_struct = nucleus_struct(qcFilter);

%create reference arrays for fluorescence heatmaps
gID_index = unique([qc_nucleus_struct.gtypeID],'stable');
gType_index = unique({qc_nucleus_struct.genotype}, 'stable');
timeVec = nan(1,numel([qc_nucleus_struct.time_interp]));
apVec = nan(1,numel([qc_nucleus_struct.ap_vector_interp]));
fluoVec = nan(1,numel([qc_nucleus_struct.fluo_interp]));
gIDVec = nan(1,numel([qc_nucleus_struct.time_interp]));

%create reference arrays for fractional activity
fonVec = nan(1,numel(qc_nucleus_struct));
foffVec = nan(1,numel(qc_nucleus_struct));
fmeanAPVec = nan(1,numel(qc_nucleus_struct));
fgIDVec = nan(1,numel(qc_nucleus_struct));

first = 1;
for n = 1:numel(qc_nucleus_struct)
    last = first+numel(qc_nucleus_struct(n).time_interp)-1;
    timeVec(first:last) = (qc_nucleus_struct(n).time_interp)/60; %want in terms of minutes
    apVec(first:last) = (qc_nucleus_struct(n).ap_vector_interp)*100;
    fluoVec(first:last) = qc_nucleus_struct(n).fluo_interp;
    gIDVec(first:last) = qc_nucleus_struct(n).gtypeID_interp;
    qcVec(first:last) = qc_nucleus_struct(n).qc;
    if qc_nucleus_struct(n).nc_trace_on_flag && qc_nucleus_struct(n).nc_trace_off_flag
       fonVec(n) = timeVec(find(~isnan([qc_nucleus_struct(n).fluo_interp]),1));
       foffVec(n) = timeVec(find(~isnan([qc_nucleus_struct(n).fluo_interp]),1, 'last')); 
    end
    fmeanAPVec(n) = qc_nucleus_struct(n).apMean*100;
    fgIDVec(n) = mean(qc_nucleus_struct(n).gtypeID);
    first = last+1;
end

%create genotype filters

wtFilter = gIDVec == gID_index(find(strcmpi(gType_index, 'Wt')));
hbFilter = gIDVec == gID_index(find(strcmpi(gType_index, 'Hb')));
gtFilter = gIDVec == gID_index(find(strcmpi(gType_index, 'Gt')));

fwtFilter = fgIDVec == gID_index(find(strcmpi(gType_index, 'Wt')));
fhbFilter = fgIDVec == gID_index(find(strcmpi(gType_index, 'Hb')));
fgtFilter = fgIDVec == gID_index(find(strcmpi(gType_index, 'Gt')));


%create arrays
%%
wtspot_array = fill_array(timeRef, apRef, wtFilter, apVec, fluoVec, 'spot', timeVec); %average fluorescence of active nuclei
wtnuc_array = fill_array(timeRef, apRef, wtFilter, apVec, fluoVec, 'nuc', timeVec); %average fluorescence of all nuclei
%%
wtfract_array = fill_array(timeRef, apRef, fwtFilter, fmeanAPVec, fluoVec, 'fract', fonVec, foffVec); %fraction of traces that start before time t and end after time t
%%
hbspot_array = fill_array(timeRef, apRef, hbFilter, apVec, fluoVec, 'spot', timeVec);
hbnuc_array = fill_array(timeRef, apRef, hbFilter, apVec, fluoVec, 'nuc', timeVec);
hbfract_array = fill_array(timeRef, apRef, fhbFilter, fmeanAPVec, fluoVec, 'fract', fonVec, foffVec);
gtspot_array = fill_array(timeRef, apRef, gtFilter, apVec, fluoVec, 'spot', timeVec);
gtnuc_array = fill_array(timeRef, apRef, gtFilter, apVec, fluoVec, 'nuc', timeVec);
%%
gtfract_array = fill_array(timeRef, apRef, fgtFilter, fmeanAPVec, fluoVec, 'fract', fonVec, foffVec);

%%
%Create individual genotype heatmaps 
spot_heatmap(wtspot_array, 'Wt', figPath)
spot_heatmap(hbspot_array, 'Hb', figPath)
spot_heatmap(gtspot_array, 'Gt', figPath)

nuc_heatmap(wtnuc_array, 'Wt', figPath)
nuc_heatmap(hbnuc_array, 'Hb', figPath)
nuc_heatmap(gtnuc_array, 'Gt', figPath)

fract_heatmap(wtfract_array, 'Wt', figPath)
fract_heatmap(hbfract_array, 'Hb', figPath)
fract_heatmap(gtfract_array, 'Gt', figPath)

%create difference arrays for wt - the two other genotypes
matrix_diff(wtspot_array, hbspot_array, 'Wt-Hb', figPath, 'spot');
matrix_diff(wtnuc_array, hbnuc_array, 'Wt-Hb', figPath, 'nuc');
matrix_diff(wtfract_array, hbfract_array,'Wt-Hb', figPath, 'fract');

matrix_diff(wtspot_array, gtspot_array, 'Wt-Gt', figPath, 'spot');
matrix_diff(wtnuc_array, gtnuc_array, 'Wt-Gt', figPath, 'nuc');
matrix_diff(wtfract_array, gtfract_array, 'Wt-Gt', figPath, 'fract');

close all

function filled_array = fill_array(timeRef, apRef, gFilter, apVec, fluoVec, varargin)
heatmap_array = NaN(numel(timeRef),numel(apRef));
spot = false;
nuc = false;
fraction = false;
timeVec = [];
onVec = [];
offVec = [];

for i = 1:length(varargin)
    if strcmpi(varargin{i},'spot')
        spot = true;
        timeVec = varargin{i+1};
    elseif strcmpi(varargin{i}, 'nuc')
        nuc = true;
        timeVec = varargin{i+1};
    elseif strcmpi(varargin{i}, 'fract')
        fraction = true;
        onVec = varargin{i+1};
        offVec = varargin{i+2};
    end
end
for a = 1:numel(apRef)
    AP = apRef(a);
    apFilter = round(apVec)==AP;
    for t = 1:numel(timeRef)
        time = timeRef(t);
        if fraction
            OnFilter = onVec < time;
            OffFilter = offVec > time;
            fractOn = sum(gFilter&OnFilter&OffFilter&apFilter)/sum(gFilter&apFilter);
            heatmap_array(t,a) = fractOn;
        elseif nuc
            tFilter = round(timeVec) == time;
            fluo = fluoVec(gFilter&apFilter&tFilter);
            fluo(isnan(fluo)) = 0;
            heatmap_array(t,a) = mean(fluo);
        elseif spot
            tFilter = round(timeVec) == time;
            fluo = fluoVec(gFilter&apFilter&tFilter);
            heatmap_array(t,a) = nanmean(fluo);
        end
    end 
end
filled_array = heatmap_array;
end

function spot_heatmap(spot_array,gName,figPath)
spot_heatmap_fig = figure;
imagesc(spot_array)
xlabel('AP position (% embryo length)')
ylabel('Time (min)')
c1 = colorbar('southoutside');
caxis([0,18e+5]); %make all the colorbars the same
c1.Label.String = 'Average Spot Intensity (au)';
title([gName, ' Spot Intensity Profile'])
saveas (spot_heatmap_fig,[figPath, 'qc_',gName, '_spot_heatmap_fig.png']);
end

function nuc_heatmap(nuc_array,gName,figPath)
nuc_heatmap_fig = figure;
imagesc(nuc_array)
xlabel('AP position (% embryo length)')
ylabel('Time (min)')
c2 = colorbar('southoutside');
caxis([0,8e+5])
%caxis to make all the colorbars the same
c2.Label.String = 'Average Nuclear Intensity (au)';
title([gName, '  Nuclear Intensity Profile'])
saveas (nuc_heatmap_fig,[figPath, 'qc_',gName, '_nuc_heatmap_fig.png']);
end

function fract_heatmap(fract_array, gName, figPath)
fract_heatmap_fig = figure;
imagesc(fract_array)
caxis([0,0.5])
xlabel('AP position (% embryo length)')
ylabel('Time (min)')
c3 = colorbar('southoutside');
c3.Label.String = 'Fraction of Active Nuclei';
title([gName, ' Fraction of Active Nuclei'])
saveas (fract_heatmap_fig, [figPath, 'qc_', gName, '_fract_heatmap_fig.png']);
end

%Create Difference Heatmaps
function matrix_diff(input1, input2, gName, figPath, varargin)
if size(input1)~=size(input2)
    error 'arrays are not the same size'
else
    for x = 1:size(input1,1)
       for y = 1:size(input1,2)
           diff_array(x,y) = input1(x,y) - input2(x,y);
       end
    end
    if strcmpi(varargin, 'spot')
        type = 'Spot Intensity Profile';
        label = 'Average Spot Intensity (au)';
        fname = ['qc_' gName '_diff_spot_fig.png'];
    elseif strcmpi(varargin, 'nuc')
        type = 'Nucleus Intensity Profile';
        label = 'Average Nuclear Intensity (au)';
        fname = ['qc_' gName '_diff_nuc_fig.png'];
    elseif strcmpi(varargin, 'fract')
        type = 'Fraction of Active Nuclei';
        label = 'Fraction of Active Nuclei';
        fname = ['qc_' gName '_diff_fract_fig.png'];
    end
    diff_heatmap_fig = figure;
    imagesc(diff_array)
    xlabel('AP position (% embryo length)')
    ylabel('Time (min)')
    c4 = colorbar('southoutside');
    c4.Label.String = label;
    title([gName,' ', type])
    saveas (diff_heatmap_fig, [figPath, fname]);
end
end