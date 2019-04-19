% Goals:
% Define a data quality flag variable that can be used to select for traces meeting a set of criteria including:
% minimum of 10 non-missing time steps DONE
% 50% non-missing observations between first and last time non-missing time
% point  DONE
% Create interpolated trace and time variables
% Interpolate everything to fall on single time grid: interpGrid = 0:20:3600
% Implement rule where all missing time points between first and last 
% observation of a trace that are more than 60 seconds from non-missing
% point are set to 0 DONE
% Next interpolate to fill in all other missing time points DONE

%run from inside src folder
clear 
close all
project = 'mse_comparison_lateralML';
currentFolder = pwd;
%specify folder for the data
slashes = strfind(currentFolder,'\');
fName = currentFolder(slashes(end)+1:end);
dataPath = ['../../dat/' project '/'];
%load the data!
load([dataPath 'nucleus_struct.mat']);

for n = 1:numel(nucleus_struct)
    %filter traces that have the desired number of non-missing time points
    minWindow = 10;
    minFill = .50; %fraction of entries between the first and last spot which must be non-missing
    nucleus_struct(n).qc = 0; %creates new quality control variable. Will be set to 1 if the trace meets the requirements
    if sum(~isnan(nucleus_struct(n).fluo)) > 0
        first_on = find(~isnan(nucleus_struct(n).fluo),1); %find the first entry that actually records a spot
        last_on = find(~isnan(nucleus_struct(n).fluo),1,'last'); %find the last entry that actially records a spot
        if sum(~isnan(nucleus_struct(n).fluo))>= minWindow %if the trace exists for 20 frames or more
            fill = sum(~isnan(nucleus_struct(n).fluo(first_on:last_on)))/(last_on - first_on); %calculates % coverage
            if fill >= minFill %if the trace has .75 non-missing entries between first and last
                nucleus_struct(n).qc = 1;
            end
        end
        %Start the interpolation
        %first, all missing time points between first and last observation
        %of a trace that are more than 60 seconds from non-missing point 
        %are set to 0
        for q = first_on:last_on
        %look forward and back 1 min to see if the point is > 60 sec away
        %from another point
            windowFilter = nucleus_struct(n).time > (nucleus_struct(n).time(q)-60) & nucleus_struct(n).time < (nucleus_struct(n).time(q)+60);
            if isnan(nucleus_struct(n).fluo(q)) && nansum(nucleus_struct(n).fluo(windowFilter)) == 0
                nucleus_struct(n).fluo(q) = 0;
            end
        end
        %linear interpolation of the NaN values between values and zeros
        nanFilter = isnan(nucleus_struct(n).fluo);
        %don't want to interpolate values before/after first/last fluorescence
        nanFilter(1:first_on-1) = 0;
        nanFilter(last_on+1:end) = 0;
        nucleus_struct(n).fluo(nanFilter) = interp1(nucleus_struct(n).time(~nanFilter), nucleus_struct(n).fluo(~nanFilter), nucleus_struct(n).time(nanFilter));
    end
    
    %now to get everything on the same time scale
    interpGrid = 0:20:3600;
    % NL: changed this so time_interp is of same length as other
    % interpolated variables
    nucleus_struct(n).interpGrid = interpGrid;
    time = nucleus_struct(n).time;
    time_interp = interpGrid(interpGrid>=time(1)&interpGrid<=time(end));
    nucleus_struct(n).time_interp = time_interp;
    interp_var_cell = {'xPosParticle', 'yPosParticle','brightestZs','fluo', 'xPos','yPos', 'ap_vector'};
    for i = 1:numel(interp_var_cell)
        varRaw = interp_var_cell{i};
        varInterp = [varRaw '_interp'];
        if numel(nucleus_struct(n).(varRaw)) >= 2
            nucleus_struct(n).(varInterp) = interp1(nucleus_struct(n).time, nucleus_struct(n).(varRaw),nucleus_struct(n).time_interp);
        else
            nucleus_struct(n).(varInterp) = []; %returns empty if there aren't enough points to interpolate
        end
    end
    
    
end
save([dataPath 'qc_nucleus_struct.mat'],'nucleus_struct')
