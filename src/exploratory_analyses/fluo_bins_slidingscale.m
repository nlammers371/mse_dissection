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
load([dataPath 'qc_nucleus_struct.mat']);

%apply quality control filtering
%qcFilter = [nucleus_struct.qc] == 1;
nucleusFilter = [nucleus_struct.nc_qc_flag] ==1;
qc_nucleus_struct = nucleus_struct(nucleusFilter);

%create reference vectors and fill them
timeVec = [];
apVec = [];
fluoVec = [];
gtypeVec = [];
setIDVec = [];
gID_index = unique([qc_nucleus_struct.gtypeID],'stable');
gType_index = unique({qc_nucleus_struct.genotype}, 'stable');
setID_index = unique([qc_nucleus_struct.setID],'stable');
apRefVec = 0:100;
timeRefVec = 0:60;

for n = 1:numel(qc_nucleus_struct)
    timeVec = [timeVec (qc_nucleus_struct(n).time_interp)/60]; %want in terms of minutes
    apVec = [apVec (qc_nucleus_struct(n).ap_vector_interp)*100];
    fluoVec = [fluoVec qc_nucleus_struct(n).fluo_interp];
    gtypeVec(numel(gtypeVec)+1:numel(gtypeVec)+numel(qc_nucleus_struct(n).time_interp)) = mean(qc_nucleus_struct(n).gtypeID);
    setIDVec(numel(setIDVec)+1:numel(setIDVec)+numel(qc_nucleus_struct(n).time_interp)) = qc_nucleus_struct(n).setID;
end
%%
%This creates an average/summation of all points over an AP bin of 3% or
%time bin of 5 min
APFigPath = [figPath 'APBins/'];
mkdir(APFigPath)
for a = 1:numel(apRefVec)-2
    apBinArrayS = NaN(numel(timeRefVec),numel(gID_index));
    apBinArrayN = NaN(numel(timeRefVec),numel(gID_index));
    apBinArrayT = NaN(numel(timeRefVec),numel(gID_index));
    apBinSEMS = NaN(numel(timeRefVec),numel(gID_index));
    apBinSEMN = NaN(numel(timeRefVec),numel(gID_index));
    apBinSEMT = NaN(numel(timeRefVec),numel(gID_index));
    
    apBinArrayF = NaN(numel(timeRefVec),numel(gID_index));
    apBinSEMF = NaN(numel(timeRefVec),numel(gID_index));
    
    apBinFilter = apVec <= (apRefVec(a)+2) & apVec > apRefVec(a);
    for t2 = 1:numel(timeRefVec)
        tFilter = round(timeVec) == timeRefVec(t2);
        for g = 1:numel(gID_index)
            gFilter = gtypeVec == gID_index(g);
            fluo = fluoVec(apBinFilter&tFilter&gFilter);
            if ~isempty(fluo)
                gsetID_index = unique(setIDVec(apBinFilter&tFilter&gFilter),'stable');
                apBinArrayS(t2,g) = nanmean(fluo);
                apBinSEMS(t2,g) = nanstd(fluo)/sqrt(sum(~isnan(fluo)));
                fluo(isnan(fluo)) = 0;
                apBinArrayN(t2,g) = mean(fluo);
                apBinSEMN(t2,g) = std(fluo)/sqrt(length(fluo));
                fluoTot = NaN(1, numel(gsetID_index));
                fluoFract = NaN(1, numel(gsetID_index));
                for s=1:numel(gsetID_index)
                    set=gsetID_index(s);
                    sFilter = setIDVec == set;
                    setFluo = fluoVec(apBinFilter&tFilter&gFilter&sFilter);
                    if ~isempty(setFluo)
                        fluoTot(s) = nansum(setFluo);
                        onFilter = setFluo > 0;
                        fluoFract(s) = sum(onFilter)/numel(setFluo);
                    else
                        fluoTot(s) = NaN;
                    end
                end
                apBinArrayT(t2,g) = nanmean(fluoTot);
                apBinSEMT(t2,g) = nanstd(fluoTot)/sqrt(sum(~isnan(fluoTot)));
                
                apBinArrayF(t2,g) = nanmean(fluoFract);
                apBinSEMF(t2,g) = nanstd(fluoFract)/sqrt(sum(~isnan(fluoFract)));
            end
        end

    end
    if sum(~isnan(apBinArrayN), 'all') ~= 0
        binName = [num2str(apRefVec(a)) '-' num2str(apRefVec(a)+2)];
        nmean_fluo_fig = figure('visible','off');
        for k = 1:size(apBinArrayN,2)
            errorbar(timeRefVec,apBinArrayN(:,k),apBinSEMN(:,k),'Marker','.', 'LineWidth', 1.5)
            hold on
        end
        hold off
        % make labels
        xlabel('Time(min)')
        ylabel('Average Fluorescence (au)')
        title(['Average nuclear fluorescence, AP bin: ' binName '% embryo length'])
        ylim([0,8e5])
        xlim([0,60])
        % add grid lines
        grid on
        % make a legend
        legend(gType_index{:})
        saveas(nmean_fluo_fig,[APFigPath  'nuc_fluo_bin' binName '.png']);
        
        smean_fluo_fig = figure('visible','off');
        for k = 1:size(apBinArrayS,2)
            errorbar(timeRefVec,apBinArrayS(:,k),apBinSEMS(:,k),'Marker','.', 'LineWidth', 1.5)
            hold on
        end
        % make labels
        xlabel('Time(min)')
        ylabel('Average Fluorescence (au)')
        title(['Average spot fluorescence, AP bin: ' binName '% embryo length'])
        ylim([0,20e5])
        xlim([0,60])
        % add grid lines
        grid on
        % make a legend
        legend(gType_index{:})
        saveas(smean_fluo_fig,[APFigPath  'spot_fluo_bin' binName '.png']);
        
        tot_fluo_fig = figure('visible','off');
        for k = 1:size(apBinArrayT,2)
            errorbar(timeRefVec,apBinArrayT(:,k),apBinSEMT(:,k),'Marker','.', 'LineWidth', 1.5)
            hold on
        end
        hold off
        % make labels
        xlabel('Time(min)')
        ylabel('Total Fluorescence (au)')
        title(['Total nuclear fluorescence, AP bin: ' binName '% embryo length'])
        xlim([0,60])
        ylim([0,8e7])
        % add grid lines
        grid on
        % make a legend
        legend(gType_index{:})
        saveas(tot_fluo_fig,[APFigPath  'tot_fluo_bin' binName '.png']);
        
        fract_fluo_fig = figure('visible','off');
        for k = 1:size(apBinArrayT,2)
            errorbar(timeRefVec,apBinArrayF(:,k),apBinSEMF(:,k),'Marker','.', 'LineWidth', 1.5)
            hold on
        end
        hold off
        % make labels
        xlabel('Time(min)')
        ylabel('Fraction Active Nuclei')
        title(['Fraction of Active Nuclei, AP bin: ' binName '% embryo length'])
        xlim([0,60])
        ylim([0,0.7])
        % add grid lines
        grid on
        % make a legend
        legend(gType_index{:})
        saveas(fract_fluo_fig,[APFigPath  'fract_fluo_bin' binName '.png']);
        close all
    end
end

%%
TFigPath = [figPath '/TimeBins/'];
mkdir(TFigPath)

for t = 1:numel(timeRefVec)-5
    tBinArrayS = NaN(numel(apRefVec),numel(gID_index));
    tBinArrayN = NaN(numel(apRefVec),numel(gID_index));
    tBinArrayT = NaN(numel(apRefVec),numel(gID_index));
    tBinSEMS = NaN(numel(apRefVec),numel(gID_index));
    tBinSEMN = NaN(numel(apRefVec),numel(gID_index));
    tBinSEMT = NaN(numel(apRefVec),numel(gID_index));
    
    tBinArrayF = NaN(numel(apRefVec),numel(gID_index));
    tBinSEMF = NaN(numel(apRefVec),numel(gID_index));
    
    tBinFilter = timeVec <= timeRefVec(t)+5 & timeVec > timeRefVec(t);
    for a2 = 1:numel(apRefVec)
        apFilter = apVec > apRefVec(a2)-0.75 & apVec < apRefVec(a2) + 0.75;
        for g = 1:numel(gID_index)
            gFilter = gtypeVec == gID_index(g);
            fluo = fluoVec(tBinFilter&apFilter&gFilter);
            if ~isempty(fluo)
                gsetID_index = unique(setIDVec(apFilter&tBinFilter&gFilter),'stable');
                tBinArrayS(a2,g) = nanmean(fluo);
                tBinSEMS(a2,g) = nanstd(fluo)/sqrt(sum(~isnan(fluo)));
                fluo(isnan(fluo)) = 0;
                tBinArrayN(a2,g) = mean(fluo);
                tBinSEMN(a2,g) = std(fluo)/sqrt(length(fluo));
                fluoTot = NaN(1,numel(gsetID_index));
                fluoFract = NaN(1,numel(gsetID_index));
                for s=1:numel(gsetID_index)
                    set=gsetID_index(s);
                    sFilter = setIDVec == set;
                    setFluo = fluoVec(apFilter&tBinFilter&gFilter&sFilter);
                    if ~isempty(setFluo)
                        fluoTot(s) = nansum(setFluo);
                        onFilter = setFluo > 0;
                        fluoFract(s) = sum(onFilter)/numel(setFluo);
                    else
                        fluoTot(s) = NaN;
                    end
                end
                tBinArrayT(a2,g) = nanmean(fluoTot);
                tBinSEMT(a2,g) = nanstd(fluoTot)/sqrt(sum(~isnan(fluoTot)));
                tBinArrayF(a2,g) = nanmean(fluoFract);
                tBinSEMF(a2,g) = nanstd(fluoFract)/sqrt(sum(~isnan(fluoFract)));
            end
        end
            
    end
    if sum(~isnan(tBinArrayN), 'all') ~= 0
        binName = [num2str(timeRefVec(t)) '-' num2str(timeRefVec(t) +5)];
        nmean_fluo_fig = figure('visible','off');
        for k = 1:size(tBinArrayN,2)
            errorbar(apRefVec,tBinArrayN(:,k),tBinSEMN(:,k),'Marker','.', 'LineWidth', 1.5)
            hold on
        end
        % make labels
        xlabel('AP Position(%embryo length)')
        ylabel('Average Fluorescence (au)')
        ylim([0,8e5])
        xlim([20,55])
        title(['Average nuclear fluorescence, Time bin: ' binName 'minutes'])
        % make a legend
        legend(gType_index{:})
        % add grid lines
        grid on
        saveas(nmean_fluo_fig,[TFigPath  'nuc_fluo_bin' binName '.png']);
        
        smean_fluo_fig = figure('visible','off');
        for k = 1:size(tBinArrayS,2)
            errorbar(apRefVec,tBinArrayS(:,k),tBinSEMS(:,k),'Marker','.', 'LineWidth', 1.5)
            hold on
        end        % make labels
        xlabel('AP Position(%embryo length)')
        ylabel('Average Fluorescence (au)')
        ylim([0,20e5])
        xlim([20,55])
        title(['Average Spot fluorescence, Time bin: ' binName 'minutes'])
        % make a legend
        legend(gType_index{:})
        % add grid lines
        grid on
        saveas(smean_fluo_fig,[TFigPath  'spot_fluo_bin' binName '.png']);
        
        tot_fluo_fig = figure('visible','off');
        for k = 1:size(tBinArrayS,2)
            errorbar(apRefVec,tBinArrayT(:,k),tBinSEMT(:,k),'Marker','.', 'LineWidth', 1.5)
            hold on
        end        % make labels
        xlabel('AP Position(%embryo length)')
        ylabel('Total Fluorescence (au)')
        xlim([20,55])
        ylim([0,18e8])
        title(['Total fluorescence, Time bin: ' binName 'minutes'])
        % make a legend
        legend(gType_index{:})
        % add grid lines
        grid on
        saveas(tot_fluo_fig,[TFigPath  'tot_fluo_bin' binName '.png']);
        close all
        
        fract_fluo_fig = figure('visible','off');
        for k = 1:size(tBinArrayS,2)
            errorbar(apRefVec,tBinArrayF(:,k),tBinSEMF(:,k),'Marker','.', 'LineWidth', 1.5)
            hold on
        end        % make labels
        xlabel('AP Position(%embryo length)')
        ylabel('Fraction Active Nuclei')
        xlim([20,55])
        ylim([0,0.7])
        title(['Fraction Active Nuclei, Time bin: ' binName 'minutes'])
        % make a legend
        legend(gType_index{:})
        % add grid lines
        grid on
        saveas(fract_fluo_fig,[TFigPath  'fract_fluo_bin' binName '.png']);
        close all

    end
end
