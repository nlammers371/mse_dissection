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

%apply quality control filtering
nucleusFilter = [nucleus_struct.nc_qc_flag] == 1;
qc_nucleus_struct = nucleus_struct(nucleusFilter);

%create reference vectors and fill them
timeVec = [];
apVec = [];
fluoVec = [];
gtypeVec = [];
GiantVec = [];
KruppelVec = [];
HunchbackVec = [];
BicoidVec = [];
gID_index = unique([qc_nucleus_struct.gtypeID],'stable');
gType_index = unique({qc_nucleus_struct.genotype}, 'stable');
TF_index = {'Gt', 'Hb', 'Kr', 'Bcd'}; 
apRefVec = 0:100;
timeRefVec = 0:60;

for n = 1:numel(qc_nucleus_struct)
    timeVec = [timeVec (qc_nucleus_struct(n).time_interp)/60]; %want in terms of minutes
    apVec = [apVec (qc_nucleus_struct(n).ap_vector_interp)*100];
    fluoVec = [fluoVec qc_nucleus_struct(n).fluo_interp];
    GiantVec = [GiantVec qc_nucleus_struct(n).tf_array(:,1).'];
    KruppelVec = [KruppelVec qc_nucleus_struct(n).tf_array(:,2).'];
    HunchbackVec = [HunchbackVec qc_nucleus_struct(n).tf_array(:,3).'];
    BicoidVec = [BicoidVec qc_nucleus_struct(n).tf_array(:,4).'];
    gtypeVec(numel(gtypeVec)+1:numel(gtypeVec)+numel(qc_nucleus_struct(n).time_interp)) = mean(qc_nucleus_struct(n).gtypeID);
end

%make folders for the genotypes
for d = 1:numel(gType_index)
    gName = gType_index{d};
    APFigPath = [figPath 'APBinsTF/' gName '/'];
    mkdir(APFigPath)
    TFigPath = [figPath 'TBinsTF/' gName '/'];
    mkdir(TFigPath)
end    
%%
%Start with ap bins
for a = 1:numel(apRefVec)-2
    apBinArrayS = NaN(numel(timeRefVec),numel(gID_index));
    apBinArrayN = NaN(numel(timeRefVec),numel(gID_index));
    apBinSEMS = NaN(numel(timeRefVec),numel(gID_index));
    apBinSEMN = NaN(numel(timeRefVec),numel(gID_index));
    apTFArray = NaN(numel(timeRefVec),4);
    
    apBinFilter = apVec <= (apRefVec(a)+2) & apVec > apRefVec(a);
    for t2 = 1:numel(timeRefVec)
        tFilter = round(timeVec) == timeRefVec(t2);
        for g = 1:numel(gID_index)
            gFilter = gtypeVec == gID_index(g);
            fluo = fluoVec(apBinFilter&tFilter&gFilter);
            if ~isempty(fluo)
                apBinArrayS(t2,g) = nanmean(fluo);
                apBinSEMS(t2,g) = nanstd(fluo)/sqrt(sum(~isnan(fluo)));
                fluo(isnan(fluo)) = 0;
                apBinArrayN(t2,g) = mean(fluo);
                apBinSEMN(t2,g) = std(fluo)/sqrt(length(fluo));
            end
            Giant = GiantVec(apBinFilter&tFilter);
            Kruppel = KruppelVec(apBinFilter&tFilter);
            Hunchback = HunchbackVec(apBinFilter&tFilter);
            Bicoid = BicoidVec(apBinFilter&tFilter);
            apTFArray(t2,1) = nanmean(Giant);
            apTFArray(t2,2) = nanmean(Kruppel);
            apTFArray(t2,3) = nanmean(Hunchback);
            apTFArray(t2,4) = nanmean(Bicoid);
        end

    end
    if sum(~isnan(apBinArrayN), 'all') ~= 0
        binName = [num2str(apRefVec(a)) '-' num2str(apRefVec(a)+2)];
        %make nuclear average plots
        make_bin_plots(timeRefVec, apBinArrayN, apBinSEMN, apTFArray, binName, 'AP', 'Nuc', TF_index, gType_index,figPath)
        %make spot fluoresence plots
        make_bin_plots(timeRefVec, apBinArrayS, apBinSEMS, apTFArray, binName, 'AP', 'Spot', TF_index, gType_index, figPath)
    end
end

%%
%continue with time bins
TFigPath = [figPath '/TimeBins/'];
mkdir(TFigPath)

for t = 1:numel(timeRefVec)-5
    tTFArray = NaN(numel(apRefVec),4);
    tBinArrayS = NaN(numel(apRefVec),numel(gID_index));
    tBinArrayN = NaN(numel(apRefVec),numel(gID_index));
    tBinSEMS = NaN(numel(apRefVec),numel(gID_index));
    tBinSEMN = NaN(numel(apRefVec),numel(gID_index));
    
    tBinFilter = timeVec <= timeRefVec(t)+5 & timeVec > timeRefVec(t);
    for a2 = 1:numel(apRefVec)
        apFilter = round(apVec) == apRefVec(a2);
        for g = 1:numel(gID_index)
            gFilter = gtypeVec == gID_index(g);
            fluo = fluoVec(tBinFilter&apFilter&gFilter);
            if ~isempty(fluo)
                tBinArrayS(a2,g) = nanmean(fluo);
                tBinSEMS(a2,g) = nanstd(fluo)/sqrt(sum(~isnan(fluo)));
                fluo(isnan(fluo)) = 0;
                tBinArrayN(a2,g) = mean(fluo);
                tBinSEMN(a2,g) = std(fluo)/sqrt(length(fluo));
            end
            Giant = GiantVec(tBinFilter&apFilter);
            Kruppel = KruppelVec(tBinFilter&apFilter);
            Hunchback = HunchbackVec(tBinFilter&apFilter);
            Bicoid = BicoidVec(tBinFilter&apFilter);
            tTFArray(a2,1) = nanmean(Giant);
            tTFArray(a2,2) = nanmean(Kruppel);
            tTFArray(a2,3) = nanmean(Hunchback);
            tTFArray(a2,4) = nanmean(Bicoid);
        end
            
    end
    if sum(~isnan(tBinArrayN), 'all') ~= 0
        binName = [num2str(timeRefVec(t)) '-' num2str(timeRefVec(t)+5)];
        %make nuclear average plots
        make_bin_plots(apRefVec, tBinArrayN, tBinSEMN, tTFArray, binName, 'Time', 'Nuc', TF_index, gType_index, figPath);
        %make spot average plots
        make_bin_plots(apRefVec, tBinArrayS, tBinSEMS, tTFArray, binName, 'Time', 'Spot', TF_index, gType_index, figPath);
    end
end

function make_bin_plots(ref_vec, fluo_array, SEM_array, tf_array, binName, bin_type, spot_or_nuc, TF_index, gType_index,figPath)
%specify color order for plots
co = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0         0         0];

for k = 1:size(fluo_array,2)
    bin_fluo_fig = figure('visible','off');
    set(bin_fluo_fig,'defaultAxesColorOrder',co)
    gName = gType_index{k};
    if strcmpi(bin_type,'AP') == true
        folderName = ['APBinsTF/' gName '/'];
        xString = 'Time (min)';
        titleEnd = '% Embryo Length';
        
    elseif strcmpi(bin_type, 'Time') == true
        folderName = ['TBinsTF/' gName '/'];
        xString = 'AP Position (%embryo length)';
        titleEnd = ' minutes';
    else
        error('invalid type')
    end
    newFigPath = [figPath folderName];
    for f = 1:4
        plot(ref_vec, tf_array(:,f), 'Marker','.')
        hold on
    end
    ylabel('Trancription Factor Concentration')
    ylim([0,1.2])
    yyaxis right
    errorbar(ref_vec,fluo_array(:,k),SEM_array(:,k),'Marker','.', 'LineWidth', 1.5, 'Color', 'k')
    % make labels
    xlabel(xString)
    ylabel('Average Fluorescence (au)')
    title(['Average ' spot_or_nuc ' fluorescence, ' bin_type ' bin: ' binName titleEnd])
    if strcmpi(spot_or_nuc,'Nuc')==true 
        ylim([0,9e5])
    else
        ylim([0,20e5])
    end
    if strcmpi(bin_type,'AP') == true
        xlim([0,60])
    elseif strcmpi(bin_type, 'Time') == true
        xlim([20,60])
    end
    hold off
    % add grid lines
    grid on
    % make a legend
    legend(TF_index{:})
    saveas(bin_fluo_fig,[newFigPath  spot_or_nuc '_fluo_bin' binName '.png']);
    close all
end
end