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

%apply quality control filtering
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
apRefVec = 0:2.5:100;
maxTime = 40; %picked a time "late" in NC 14, since the paper was non-specific

for n = 1:numel(qc_nucleus_struct)
    timeVec = [timeVec (qc_nucleus_struct(n).time_interp)/60]; %want in terms of minutes
    apVec = [apVec (qc_nucleus_struct(n).ap_vector_interp)*100];
    fluoVec = [fluoVec qc_nucleus_struct(n).fluo_interp];
    gtypeVec(numel(gtypeVec)+1:numel(gtypeVec)+numel(qc_nucleus_struct(n).time_interp)) = mean(qc_nucleus_struct(n).gtypeID);
    setIDVec(numel(setIDVec)+1:numel(setIDVec)+numel(qc_nucleus_struct(n).time_interp)) = qc_nucleus_struct(n).setID;
end

timeCap = timeVec < maxTime;
%%
for g = 1:numel(gID_index)
    integVec = nan(1, numel(apRefVec));
    gFilter = gtypeVec==gID_index(g);
    for a = 1:numel(apRefVec)
        apFilter = apVec > apRefVec(a) & apVec < apRefVec(a) + 2.5;
        %at each AP position find the integrated fluorescence per set
        fluo = nan(1, numel(setID_index));
        for s = 1:numel(setID_index)
            setFilter = setIDVec == setID_index(s);
            if ~isempty(fluoVec(gFilter&apFilter&timeCap&setFilter))
                fluo(s) = nansum(fluoVec(gFilter&apFilter&timeCap&setFilter));
            end
        end
        
        if ~isempty(fluo)
            %take the average across the sets
            integVec(a) = nansum(fluo)/sum(~isnan(fluo));
        end
    end
    embryo_fig = figure(g);
    heatmap(integVec,'CellLabelColor','none','MissingDataColor',[1 1 1],...
        'ColorLimits',[0,5.5e8],'GridVisible','off')
    saveas(embryo_fig, [figPath gType_index{g} '_pseudo_embryo.png']);
    
end