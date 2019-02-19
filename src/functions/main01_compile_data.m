% main01_compile_data(project, FolderPath, keyword, include_vec)
%
% DESCRIPTION
% Funcion to compile relevant outputs from image analysis pipeline across
% multiple experiments

%
% ARGUMENTS
% project: master ID variable 
%
% FolderPath: Full or relative path to DynamicsResults folder (pr
%               equivalent)
%
% keyword: String contained in all folder names for projects one wishes to 
%          compile. For instance: 'Eve2MS2', will full all projects
%          containing this string
%
% OPTIONS
% include_vec:  Vector specifying IDs within list of projects
%               pulled by 'keyword' to keep. If passed as empty vector, all
%               matching projects will be taken
% DropboxFolder: Pass this option, followed by the path to data folder 
%                where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
% first_nc: script defaults to taking only nc14 traces. If you want
%           earlier traces, pass 'first_nc', followed by desired nuclear cycle
%           number
% exp_type: by default, code assumes input/output experiment. If no
%           protein channel exists, pass 'exp_type', follwed by any
%           other string. In future can add additional
%           functionalities for other exp types
%
% OUTPUT: nucleus_struct: compiled data set contain key nucleus and
% particle attributes

function nucleus_struct = main01_compile_data(project,FolderPath,keyword,varargin)

% set defaults
includeVec = [];
firstNC = 14;
expType = 'input_output';
dataPath = ['../dat/' project '/']; % data mat directory
for i = 1:numel(varargin)
    if ischar(varargin{i})
        if ismember(varargin{i},{'includeVec','firstNC','expType'})
           eval([varargin{i} '=varargin{i+1};']);
        end
        if strcmpi(varargin{i},'dropboxFolder')
            dataPath = [dropboxFolder '\ProcessedEnrichmentData\' project '/'];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Set Path Specs, ID Vars %%%%%%%%%%%%%%%%%%%%%%%%

% folders

%%% make filepath
mkdir(dataPath);
%%% assign save names
nucleus_name = [dataPath 'nucleus_struct.mat']; % names for compiled elipse struct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Obtain Relevant Filepaths %%%%%%%%%%%%%%%%%%%%%%%
% obtain and store set names
dirinfo = dir([FolderPath '*' keyword '*']);
dirinfo(~[dirinfo.isdir]) = []; %remove non-directories
if isempty(includeVec)
    includeVec = 1:numel(dirinfo);
end
% initialize filename vectors
cp_filenames = {}; % particles
cn_filenames = {}; % protein data
ap_filenames = {}; % ap info
nc_filenames = {}; % nuclei
fov_filenames = {}; % fov info
sp_filenames = {}; % spot info
pa_filenames = {}; % particle info that contains z info
for d = 1:numel(includeVec)
    thisdir = dirinfo(includeVec(d)).name;            
    % append file paths
    cn_filenames = [cn_filenames {[thisdir '/CompiledNuclei.mat']}];
    cp_filenames = [cp_filenames {[thisdir '/CompiledParticles.mat']}];    
    ap_filenames = [ap_filenames {[thisdir '/APDetection.mat']}];    
    nc_filenames = [nc_filenames {[thisdir '/' thisdir '_lin.mat']}];           
    fov_filenames = [fov_filenames {[thisdir '/FrameInfo.mat']}];
    sp_filenames = [sp_filenames {[thisdir '/Spots.mat']}];
    pa_filenames = [pa_filenames {[thisdir '/Particles.mat']}];
end

% generate set key data structure
set_key = array2table(includeVec','VariableNames',{'setID'});
dirs = {dirinfo(includeVec).name};
set_key.prefix = dirs';
save([dataPath 'set_key.mat'],'set_key')

% Generate master structure with info on all nuclei and traces in
% constituent sets
nucleus_struct = [];
% Loop through filenames    
for i = 1:length(cp_filenames) 
    % read in raw files
    ap_flag = 1;
    try
        load([FolderPath ap_filenames{i}]) % AP Info   
    catch
        warning('No AP Data detected. Proceeding without AP info')
        ap_flag = 0;
    end
    load([FolderPath nc_filenames{i}]) % Ellipse Info
    load([FolderPath fov_filenames{i}]) % FrameInfo Info
    load([FolderPath sp_filenames{i}]) % Spots info
    load([FolderPath pa_filenames{i}]); % Raw Particles
    
    % extract data structures
    processed_data = load([FolderPath cp_filenames{i}]); % processed particles    
    input_output = 0;
    if strcmpi(expType, 'input_output')
        input_output = 1;
        protein_data = load([FolderPath cn_filenames{i}]); % processed nulcei 
    end
    
    % set identifier
    setID = includeVec(i);            
    % pull trace and nucleus variables
    time_raw = processed_data.ElapsedTime*60; % time vector            
    traces_raw = processed_data.AllTracesVector; % array with a column for each trace 
    % check to see if traces are stored in cell array
    if iscell(traces_raw)
        traces_raw = traces_raw{1};
    end
    frames_raw = 1:length(time_raw); % Frame list   
    
    first_frame = max([1, processed_data.(['nc' num2str(firstNC)])]); % Default to start of nc14 for the moment    
    % filter trace mat and time
    traces_clean = traces_raw(first_frame:end,:);
    time_clean = time_raw(first_frame:end);    
    time_clean = time_clean - min(time_clean); % Normalize to start of nc14    
    frames_clean = frames_raw(first_frame:end);    
    
    if input_output
        % extract protein data
        cp_protein = protein_data.CompiledNuclei;   
        if iscell(cp_protein)
            cp_protein= cp_protein{1};
        end
    end
    
    % compile schnitz info
    s_cells = struct;
    e_pass = 1;    
    for e = 1:length(schnitzcells)
        e_frames = schnitzcells(e).frames;
        nc_filter = ismember(e_frames,frames_clean);
        nc_frames = e_frames(nc_filter);
        if length(nc_frames) >= 1 % skip nuclei not desired nc range                     
            %Will be set to particle real values for nuclei with matching
            %particle
            s_cells(e_pass).ParticleID = NaN;
            s_cells(e_pass).xPosParticle = NaN(1,sum(nc_filter));
            s_cells(e_pass).yPosParticle = NaN(1,sum(nc_filter));
            s_cells(e_pass).brightestZs = NaN(1,sum(nc_filter));
            s_cells(e_pass).spot_z_cell = cell(1,sum(nc_filter));
            s_cells(e_pass).fluo = NaN(1,sum(nc_filter));
            
            % add core nucleus info
            x = schnitzcells(e).cenx;            
            y = schnitzcells(e).ceny;                           
            s_cells(e_pass).xPos = x(nc_filter);
            s_cells(e_pass).yPos = y(nc_filter); 
            s_cells(e_pass).frames = nc_frames';            
            s_cells(e_pass).Nucleus = e;     
            if ap_flag
                s_cells(e_pass).ap_vector = schnitzcells(e).APpos(nc_filter);
                s_cells(e_pass).apMean = mean(schnitzcells(e).APpos(nc_filter));
            end
            s_cells(e_pass).ncID = eval([num2str(setID) '.' sprintf('%04d',e)]);
            s_cells(e_pass).xMean = mean(x(nc_filter));
            s_cells(e_pass).yMean = mean(y(nc_filter));
            s_cells(e_pass).ncStart = firstNC;
            % time and set info
            s_cells(e_pass).time = time_clean(ismember(frames_clean,nc_frames));                        
            s_cells(e_pass).setID = setID;
            fn = cp_filenames{i}; % Get filename to store in struct  
            fn = fn(1:strfind(fn,'/')-1);
            s_cells(e_pass).source_path = fn;                        
            s_cells(e_pass).PixelSize = FrameInfo(1).PixelSize;       
            s_cells(e_pass).ap_flag = ap_flag;
            % add protein info      
            if input_output
                prot_nuc = cp_protein([cp_protein.schnitz] == e);                     
                pt_vec = prot_nuc.FluoMax(nc_filter);
                s_cells(e_pass).protein = pt_vec';
            end
            e_pass = e_pass + 1;
        end
    end
    
    % now add particle info
    
    % Index vector to cross-ref w/ particles            
    e_index = [s_cells.Nucleus]; 
    
    % Extract compiled particles structure       
    cp_particles = processed_data.CompiledParticles; 
    if iscell(cp_particles)
        cp_particles = cp_particles{1};
    end    
    
    % iterate through traces
    for j = 1:size(traces_clean,2)  
        % Raw fluo trace
        raw_trace = traces_clean(:,j); 
        % Get nucleus ID
        schnitz = cp_particles(j).schnitz;        
        % skip particles not in nc range
        if sum(~isnan(raw_trace)) == 0
            continue
        end
        trace_start = find(~isnan(raw_trace),1);
        trace_stop = find(~isnan(raw_trace),1,'last');        
        %Creat versions with all intervening frames present (missing frames
        %appear as NaNs)
        trace_full = raw_trace(trace_start:trace_stop)';                                   
        frames_full = frames_clean(trace_start:trace_stop);                       
        % Find intersection btw full frame range and CP frames        
        raw_pt_frames = cp_particles(j).Frame;
        cp_frames = raw_pt_frames(ismember(raw_pt_frames,frames_full));        
        
        % Initialize arrays to store Z info
        bZ = NaN(1, length(frames_full));
        z_cell = cell(1, length(frames_full));
        % Get spot Z location info
        part_id = cp_particles(j).OriginalParticle;        
        particle_frames_raw = Particles(part_id).Frame;
        for id = 1:length(cp_frames)            
            abs_idx = find(cp_frames(id) == frames_full);
            rel_id = raw_pt_frames==cp_frames(id);
            if Particles(part_id).Index(particle_frames_raw==cp_frames(id)) > length(Spots(cp_frames(id)) ...
                    .Fits)
                error('Mismatch between Particles and Spots dimensions')                                
            else
                bZ(abs_idx) = Spots(cp_frames(id)) ...
                        .Fits(Particles(part_id).Index(rel_id)).brightestZ;
                z_cell{abs_idx} = Spots(Particles(part_id).Frame(rel_id)) ...
                        .Fits(Particles(part_id).Index(rel_id)).z;
            end
        end
        % find corresponding nucleus index        
        nc_ind = find(e_index==schnitz);        
        if length(nc_ind) ~= 1
            error('Problem with Particle-Nucleus Crossref')
        end                                         
        % find overlap between nucleus and trace
        nc_frames = s_cells(nc_ind).frames;         
        spot_filter = ismember(nc_frames,frames_full);            
        s_cells(nc_ind).spot_frames = spot_filter;
        if sum(spot_filter) < numel(frames_full)
            error('Inconsistent particle and nucleus frames')
        end
        % record fluorescence info             
        s_cells(nc_ind).fluo(spot_filter) = trace_full;               
        % x and y info                                
        s_cells(nc_ind).xPosParticle(ismember(nc_frames,cp_frames)) = ...
            cp_particles(j).xPos(ismember(cp_frames,nc_frames));
        s_cells(nc_ind).yPosParticle(ismember(nc_frames,cp_frames)) = ...
            cp_particles(j).yPos(ismember(cp_frames,nc_frames));
        % add z info                       
        s_cells(nc_ind).brightestZs(spot_filter) = bZ;   
        s_cells(nc_ind).spot_z_cell(spot_filter) = z_cell;
        % Identifier variables                        
        particle = cp_particles(j).OriginalParticle;                        
        s_cells(nc_ind).ParticleID = eval([num2str(setID) '.' sprintf('%04d',particle)]);        
    end      
    nucleus_struct = [nucleus_struct  s_cells];        
end
% save
save(nucleus_name ,'nucleus_struct') 

