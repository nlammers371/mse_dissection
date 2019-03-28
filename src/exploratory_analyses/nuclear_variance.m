clear
dataPath = '/Users/annikamartin/Documents/Research/Rotation3/Exploratory_Images/dat/mse_comparison_mcp2f/';
figPath = '/Users/annikamartin/Documents/Research/Rotation3/Exploratory_Images/fig/mse_comparison_mcp2f/';
load([dataPath 'qc_nucleus_struct.mat']);

%reference vectors
timeRef = 1:60;
apRef = 1:100;

%apply quality control filtering
qcFilter = [nucleus_struct.qc] == 1;
qc_nucleus_struct = nucleus_struct(qcFilter);

%create reference arrays

for n = 1:numel(qc_nucleus_struct)
    timeVec(n,:) = qc_nucleus_struct(n).time_interp; %structure format (nucleus:time steps)
    apVec(n,:) = (qc_nucleus_struct(n).ap_vector_interp)*100;
    fluoVec(n,:) = qc_nucleus_struct(n).fluo_interp;
    gtypeVec(n,1:181) = {qc_nucleus_struct(n).genotype};
    onVec(n,1:181) = timeVec(n,find(~isnan([qc_nucleus_struct(n).fluo_interp]),1));
    offVec(n,1:181) = timeVec(n,find(~isnan([qc_nucleus_struct(n).fluo_interp]),1, 'last'));
end

%create genotype filters
wtFilter = strcmp(gtypeVec, 'Wt');
hbFilter = strcmp(gtypeVec, 'Hb');
gtFilter = strcmp(gtypeVec, 'Gt');

%create arrays
wtspot_array = NaN(numel(timeRef),numel(apRef)); %average fluorescence of active nuclei
wtnuc_array = NaN(numel(timeRef),numel(apRef)); %average fluorescence of all nuclei
wtfract_array = NaN(numel(timeRef),numel(apRef)); %fraction of traces that start before time t and end after time t
hbspot_array = NaN(numel(timeRef),numel(apRef));
hbnuc_array = NaN(numel(timeRef),numel(apRef));
hbfract_array = NaN(numel(timeRef),numel(apRef));
gtspot_array = NaN(numel(timeRef),numel(apRef));
gtnuc_array = NaN(numel(timeRef),numel(apRef));
gtfract_array = NaN(numel(timeRef),numel(apRef));

%fill arrays
for a = 1:numel(apRef)
    AP = apRef(a);
    apFilter = round(apVec(:,:))==AP;
    for t = 1:numel(timeRef)
        time = timeRef(t)*60;
        tFilter = timeVec(:,:)==time;
        wtFluo = fluoVec(wtFilter&apFilter&tFilter);
        wtspot_array(t,a) = nanmean(wtFluo);
        hbFluo = fluoVec(hbFilter&apFilter&tFilter);
        hbspot_array(t,a) = nanmean(hbFluo);
        gtFluo = fluoVec(gtFilter&apFilter&tFilter);
        gtspot_array(t,a) = nanmean(gtFluo);
        
        wtFluo(isnan(wtFluo))=0;
        wtnuc_array(t,a) = mean(wtFluo);
        hbFluo(isnan(hbFluo))=0;
        hbnuc_array(t,a) = mean(hbFluo);
        gtFluo(isnan(gtFluo))=0;
        gtnuc_array(t,a) = mean(gtFluo);
        
        onFilter = onVec(:,:)<time;
        offFilter = offVec(:,:)>time;
        wtnucOn = numel(fluoVec(wtFilter&onFilter&offFilter&apFilter&tFilter));
        wtnucTot = numel(wtFluo);
        wtfract_array(t,a) = wtnucOn/wtnucTot;
        hbnucOn = numel(fluoVec(hbFilter&onFilter&offFilter&apFilter&tFilter));
        hbnucTot = numel(hbFluo);
        hbfract_array(t,a) = hbnucOn/hbnucTot;
        gtnucOn = numel(fluoVec(gtFilter&onFilter&offFilter&apFilter&tFilter));
        gtnucTot = numel(gtFluo);
        gtfract_array(t,a) = gtnucOn/gtnucTot;
        
    end
    
end

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

function spot_heatmap(spot_array,gName,figPath)
spot_heatmap_fig = figure;
imagesc(spot_array)
colorbar
caxis([0,20e+5]); %make all the colorbars the same
xlabel('AP position (% embryo length)')
ylabel('Time (min)')
c1 = colorbar;
c1.Label.String = 'Average Spot Intensity (au)';
title([gName, ' Spot Intensity Profile'])
saveas (spot_heatmap_fig,[figPath, gName, '_spot_heatmap_fig.png']);
end

function nuc_heatmap(nuc_array,gName,figPath)
nuc_heatmap_fig = figure;
imagesc(nuc_array)
colorbar
caxis([0,6e+5])
%caxis to make all the colorbars the same
xlabel('AP position (% embryo length)')
ylabel('Time (min)')
c2 = colorbar;
c2.Label.String = 'Average Nuclear Intensity (au)';
title([gName, '  Nuclear Intensity Profile'])
saveas (nuc_heatmap_fig,[figPath, gName, '_nuc_heatmap_fig.png']);
end

function fract_heatmap(fract_array, gName, figPath)
fract_heatmap_fig = figure;
imagesc(fract_array)
colorbar
caxis([0,1])
xlabel('AP position (% embryo length)')
ylabel('Time (min)')
c3 = colorbar;
c3.Label.String = 'Fraction of Active Nuclei';
title([gName, ' Fraction of Active Nuclei'])
saveas (fract_heatmap_fig, [figPath, gName, '_fract_heatmap_fig.png']);
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
        fname = [gName '_diff_spot_fig.png'];
    elseif strcmpi(varargin, 'nuc')
        type = 'Nucleus Intensity Profile';
        label = 'Average Nuclear Intensity (au)';
        fname = [gName '_diff_nuc_fig.png'];
    elseif strcmpi(varargin, 'fract')
        type = 'Fraction of Active Nuclei';
        label = 'Fraction of Active Nuclei';
        fname = [gName '_diff_fract_fig.png'];
    end
    diff_heatmap_fig = figure;
    imagesc(diff_array)
    colorbar
    xlabel('AP position (% embryo length)')
    ylabel('Time (min)')
    c4 = colorbar;
    c4.Label.String = label;
    title([gName,' ', type])
    saveas (diff_heatmap_fig, [figPath, fname]);
end
end