%Run from inside the exploratory_analyses folder
clear
close all
project = 'mse_comparison_lateralML';
currentFolder = pwd;
%find file name and make another file of the same name in the data folder
%to hold our figures
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
figPath = ['../../fig/' project '/' fName '/OnOffTimes/'];
mkdir(figPath);
%specify where the data is coming from and load the data
dataPath = ['../../dat/' project '/'];
load([dataPath 'qc_nucleus_struct.mat']);

%apply nucleus quality control filtering
qcFilter = [nucleus_struct.qc] == 1;
qc_nucleus_struct = nucleus_struct(qcFilter);

%create empty vectors
onTimeVec = NaN(1,numel(qc_nucleus_struct)); %time of first fluo
offTimeVec = NaN(1,numel(qc_nucleus_struct)); %time of last fluo
windowVec =  NaN(1,numel(qc_nucleus_struct));
meanAPVec = NaN(1,numel(qc_nucleus_struct)); %mean AP position of nucleus
gIDVec = NaN(1,numel(qc_nucleus_struct));
setIDVec = NaN(1,numel(qc_nucleus_struct));

gID_index = unique([qc_nucleus_struct.gtypeID],'stable');
gType_index = unique({qc_nucleus_struct.genotype}, 'stable');
setID_index = unique([qc_nucleus_struct.setID], 'stable');
lastFrame_index = NaN(1,numel(setID_index));

%calculate the final time sampled for each set
for s = 1:numel(setID_index)
    setFilter = [qc_nucleus_struct.setID] == setID_index(s);
    lastFrame_index(1,s) = max([qc_nucleus_struct(setFilter).time]);
end

%fill vectors
for i=1:numel(qc_nucleus_struct)
   meanAPVec(i) = qc_nucleus_struct(i).apMean*100;
   gIDVec(i) = mean(qc_nucleus_struct(i).gtypeID);
   setIDVec(i) = qc_nucleus_struct(i).setID;
   if qc_nucleus_struct(i).nc_trace_on_flag == true
       onTime = qc_nucleus_struct(i).time(find(~isnan(qc_nucleus_struct(i).fluo),1));
   else
       onTime = NaN;
   end
   if qc_nucleus_struct(i).nc_trace_off_flag == true
       offTime = qc_nucleus_struct(i).time(find(~isnan(qc_nucleus_struct(i).fluo), 1, 'last'));
   else
       offTime = NaN;
   end
   
   onTimeVec(i) = onTime;
   offTimeVec(i) = offTime;
   windowVec(i) = (offTime-onTime)/60;
end

%now sort by AP position and genotype
APRefVec = 20:60;
avgOnVec = NaN(numel(gID_index), numel(APRefVec));
avgOffVec = NaN(numel(gID_index), numel(APRefVec));
SEMOnVec = NaN(numel(gID_index), numel(APRefVec));
SEMOffVec = NaN(numel(gID_index), numel(APRefVec));
avgWindowVec = NaN(numel(gID_index), numel(APRefVec));
SEMWindowVec = NaN(numel(gID_index), numel(APRefVec));
for g = 1:numel(gID_index)
    gFilter = gIDVec == gID_index(g);
    for a = 1:numel(APRefVec)
        apFilter = round(meanAPVec)==APRefVec(a);
        onTime = nanmean(onTimeVec(apFilter&gFilter))/60; %convert to minutes
        offTime = nanmean(offTimeVec(apFilter&gFilter))/60;
        onSEM = nanstd(onTimeVec(apFilter&gFilter)/60)/sqrt(sum(~isnan(onTimeVec(apFilter&gFilter))));
        offSEM = nanstd(offTimeVec(apFilter&gFilter)/60)/sqrt(sum(~isnan(offTimeVec(apFilter&gFilter))));
        meanWindow = nanmean(windowVec(apFilter&gFilter));
        windowSEM = nanstd(windowVec(apFilter&gFilter))/sqrt(sum(~isnan(windowVec(apFilter&gFilter))));
        avgOnVec(g,a) = onTime;
        avgOffVec(g,a) = offTime;
        SEMOnVec(g,a) = onSEM;
        SEMOffVec(g,a) = offSEM;
        avgWindowVec(g,a) = meanWindow;
        SEMWindowVec(g,a) = windowSEM;
    end
    on_time_fig = figure;
    scatter(meanAPVec(gFilter),(onTimeVec(gFilter)/60), 'MarkerEdgeAlpha', 0.5)
    hold on
    errorbar(APRefVec,avgOnVec(g,:),SEMOnVec(g,:), 'Marker', '.', 'LineWidth', 1.25)
    xlabel('AP position (%embryo length)')
    ylabel('Average Nuclear On Time (minutes)')
    ylim([10,60])
    xlim([20,50])
    title ([gType_index{g} ' Average Nuclear On Time'])
    grid on
    hold off
    saveas (on_time_fig,[figPath gType_index{g} '_on_time_fig.png']);
    
    off_time_fig = figure;
    scatter(meanAPVec(gFilter),(offTimeVec(gFilter)/60), 'MarkerEdgeAlpha', 0.5)
    hold on
    errorbar(APRefVec,avgOffVec(g,:), SEMOffVec(g,:), 'Marker', '.', 'LineWidth', 1.25)
    xlabel('AP position (%embryo length)')
    ylabel('Average Nuclear Off Time (minutes)')
    ylim([10,60])
    xlim([20,50])
    title ([gType_index{g} ' Average Nuclear Off Time'])
    grid on
    hold off
    saveas (off_time_fig,[figPath gType_index{g} '_off_time_fig.png']);
    
    window_fig = figure;
    scatter(meanAPVec(gFilter),(windowVec(gFilter)), 'MarkerEdgeAlpha', 0.5)
    hold on
    errorbar(APRefVec,avgWindowVec(g,:), SEMWindowVec(g,:), 'Marker', '.', 'LineWidth', 1.25)
    xlabel('AP position (%embryo length)')
    ylabel('Average Transcription Window (minutes)')
    ylim([0,40])
    xlim([20,50])
    title ([gType_index{g} ' Average Transcriptional Window'])
    grid on
    hold off
    saveas (window_fig,[figPath gType_index{g} '_window_fig.png']);
end
%%
on_all_fig = figure;
for i = 1:numel(gID_index)
    errorbar(APRefVec,avgOnVec(i,:), SEMOnVec(i,:), 'Marker', '.', 'LineWidth', 1.25)
    hold on
end
xlabel('AP position (%embryo length)')
ylabel('Average Nuclear On Time (minutes)')
ylim([10,60])
xlim([20,50])
legend(gType_index{:})
title ('Average Nuclear On Time')
grid on
hold off
saveas (on_all_fig,[figPath 'Combined_on_time_fig.png']);

off_all_fig = figure;
for j = 1:numel(gID_index)
    errorbar(APRefVec,avgOffVec(j,:), SEMOffVec(j,:), 'Marker', '.', 'LineWidth', 1.25)
    hold on
end
xlabel('AP position (%embryo length)')
ylabel('Average Nuclear Off Time (minutes)')
ylim([10,60])
xlim([20,50])
legend(gType_index{:})
title ('Average Nuclear Off Time')
grid on
hold off
saveas (off_all_fig,[figPath 'Combined_off_time_fig.png']);

window_all_fig = figure;
for j = 1:numel(gID_index)
    errorbar(APRefVec,avgWindowVec(j,:), SEMWindowVec(j,:), 'Marker', '.', 'LineWidth', 1.25)
    hold on
end
xlabel('AP position (%embryo length)')
ylabel('Average Transcriptional Window (minutes)')
ylim([0,40])
xlim([20,50])
legend(gType_index{:})
title ('Average Transcriptional Window')
grid on
hold off
saveas (off_all_fig,[figPath 'Combined_window_fig.png']);
