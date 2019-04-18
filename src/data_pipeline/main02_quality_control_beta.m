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

function nucleus_struct = main02_quality_control_beta(project,varargin)
%%% AM: Still need to add varargin implementation(minwindow, minfill, etc.)

dataPath = ['../../dat/' project '/'];
%load the data!
load([dataPath 'nucleus_struct.mat'],'nucleus_struct');

for n = 1:numel(nucleus_struct)
    %filter traces that have the desired number of non-missing time points
    minWindow = 10;
    minFill = .50; %fraction of entries between the first and last spot which must be non-missing
    nucleus_struct(n).qc = false; %creates new quality control variable. Will be set to 1 if the trace meets the requirements
    if sum(~isnan(nucleus_struct(n).fluo)) > 0
        first_on = find(~isnan(nucleus_struct(n).fluo),1); %find the first entry that actually records a spot
        last_on = find(~isnan(nucleus_struct(n).fluo),1,'last'); %find the last entry that actially records a spot
        if sum(~isnan(nucleus_struct(n).fluo))>= minWindow %if the trace exists for 10 frames or more
            fill = sum(~isnan(nucleus_struct(n).fluo(first_on:last_on)))/(last_on - first_on); %calculates % coverage
            if fill >= minFill %if the trace has 50% non-missing entries between first and last
                nucleus_struct(n).qc = true;
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
    interp_var_cell = {'xPosParticle', 'yPosParticle','brightestZs','fluo', 'xPos','yPos', 'ap_vector', 'gtypeID'};
    %AM: Added gtypeID variable to make genotype filtering easier
    if numel(nucleus_struct(n).time) >= 2
        for i = 1:numel(interp_var_cell)
            varRaw = interp_var_cell{i};
            varInterp = [varRaw '_interp'];        
            nucleus_struct(n).(varInterp) = interp1(nucleus_struct(n).time, nucleus_struct(n).(varRaw),nucleus_struct(n).time_interp);
        end
    else    
        [~, mi] = min(abs(interpGrid - nucleus_struct(n).time));
        nucleus_struct(n).time_interp = interpGrid(mi);
        for i = 1:numel(interp_var_cell)
            varRaw = interp_var_cell{i};
            varInterp = [varRaw '_interp'];        
            nucleus_struct(n).(varInterp) = nucleus_struct(n).(varRaw);
        end
    end  
    
    %%% NL: making nucleus-based qc flags
    %  AM: added condition that nucleus must be in frame when it turns off
    nc_qc_flag = numel(nucleus_struct(n).time_interp) > 20; 
    nc_trace_off_flag = false;
    nc_trace_on_flag = false;
    fluo_interp = nucleus_struct(n).fluo_interp;
    
    if any(~isnan(fluo_interp))
        trace_start = find(~isnan(fluo_interp),1);
        trace_stop = find(~isnan(fluo_interp),1,'last');
        last_frame = find(fluo_interp,1,'last');
        nc_trace_on_flag = nucleus_struct(n).qc &  trace_start >= 6 & nc_qc_flag;
        if trace_stop == last_frame
            setID = nucleus_struct(n).setID;
            setFilter = [nucleus_struct.setID] == setID;
            max_time = max([nucleus_struct(setFilter).time_interp]); %defines the final frame taken in this imaging set
            if fluo_interp(last_frame) == max_time
                %for now, if the spot is still visible at the end of
                %imaging, the off time counts
                nc_trace_off_flag = true; 
            end
        else
            nc_trace_off_flag = nucleus_struct(n).qc & trace_stop < last_frame & nc_qc_flag;
        end            
    end
    nucleus_struct(n).nc_qc_flag = nc_qc_flag;
    nucleus_struct(n).nc_trace_on_flag = nc_trace_on_flag;
    nucleus_struct(n).nc_trace_off_flag = nc_trace_off_flag;        
end
save([dataPath 'qc_nucleus_struct.mat'],'nucleus_struct')
