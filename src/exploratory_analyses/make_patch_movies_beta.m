function make_patch_movies_beta(Prefix, PreProcPath, PostProcPath, nc, varargin)
%
% Script to generate movies with overlaid colored patches indicating mean activity
% or EVER ON status
%
%author: Nick Lammers
%
%
%TO DO: 
% 1. Add instantaneous fluorescence visualization
% 2. Scale spot channel and opacity better 
%%% Load Data & Create Write Paths

close all;

fractionOn = 0; % if 1,  only indicate whether nucleus turns on. 0 Mean rate patches
meanRate = 1;
instantRate = 0;
visible = 0; % if 1, shows each frame
saveFigs = 1;
last_time = NaN;
set_color = 0;
for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'on_only')
        fractionOn = 1;
        meanRate = 0;
    elseif strcmpi(varargin{i}, 'instantRate')
        instantRate = 1;
        meanRate = 0;
    elseif strcmpi(varargin{i}, 'visible')
        visible = 1;
    elseif strcmpi(varargin{i}, 'noSave')
        saveFigs = 0;
    elseif strcmpi(varargin{i}, 'last_time')
        last_time = varargin{i+1};
    elseif strcmpi(varargin{i}, 'set_color')
        set_color = 1;
        fill_color = [varargin{i+1}]/256;
    end
end


% [~,~,PrefixDropboxFolder,~, PreProcPath,...
%     ~, Prefix, ~,~,~,~, ~] = readMovieDatabase(Prefix);
if ~set_color
    if fractionOn
        fill_color = [234 194 100]/256; % patch color
    elseif meanRate 
        fill_color = [85 169 116]/256; % patch color
    elseif instantRate 
        fill_color = [115 143 193]/256; % patch color
    end
end
% make write path
PostProcPath = [PostProcPath '/' Prefix '/'];
PreProcPath = [PreProcPath '/' Prefix '/'];
if fractionOn
    WritePath = [PostProcPath '/FractionOnFrames/'];
    save_prefix = 'fraction_on';
elseif meanRate
    WritePath = [PostProcPath '/MeanRateFrames/'];
    save_prefix = 'mean_rate';
elseif instantRate
    WritePath = [PostProcPath '/InstantRateFrames/'];
    save_prefix = 'instant_rate';
else
    error('no patch option selected')
end
mkdir(WritePath);

%Load the data
load([PostProcPath,'\CompiledParticles.mat']);
load([PostProcPath,'\' Prefix,'_lin.mat'], 'schnitzcells');
load([PostProcPath,'\Ellipses.mat'], 'Ellipses');
load([PostProcPath,'\FrameInfo.mat'], 'FrameInfo');

%%% calculate additional (derivative) movie parameters
PixelSize = FrameInfo(1).PixelSize;
MaxRadius = 5 / PixelSize; % um
xDim = FrameInfo(1).PixelsPerLine;
yDim = FrameInfo(1).LinesPerFrame;
zSlices = FrameInfo(1).NumberSlices;

[px, py] = meshgrid(1:xDim,1:yDim);
channel = 1; %this can be modified for multiple channels if needed

if iscell(AllTracesVector)
    AllTracesVector = AllTracesVector{channel};
end
if iscell(CompiledParticles)
    CompiledParticles = CompiledParticles{channel};
end

MaxmRNAFluo = prctile(AllTracesVector(:),99);    
first_frame = eval(['nc' num2str(nc)]);

if ~isnan(last_time)
    last_frame = find(floor(ElapsedTime-ElapsedTime(first_frame)) == last_time,1);
elseif nc == 14
    last_frame = numel(ElapsedTime);
else
    try
        last_frame = eval(['nc' num2str(nc+1)]);
    catch
        last_frame = numel(ElapsedTime);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate Movies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if visible
    OverlayFig = figure();
else
    OverlayFig = figure('Visible','off');
end

OverlayAxes = axes(OverlayFig);


FrameRange=first_frame:last_frame;

%Iterate Through Frames
for Dummy = FrameRange    
    CurrentFrame = max(FrameRange) - Dummy + first_frame;
   
    %Track pixel assignments
    NucleusStateMat = zeros(yDim,xDim,3);
    NucleusIDMat = zeros(yDim,xDim);
    NucleusDistMat = ones(yDim,xDim)*MaxRadius;
    %Loop through ALL nuclei (inactive and active) and assign patches
    for s = 1:length(schnitzcells)
        MaxFrame = max(schnitzcells(s).frames);
        MinFrame = min(schnitzcells(s).frames);        
        CurrentEllipse= schnitzcells(s).cellno(...
                        schnitzcells(s).frames==...
                        CurrentFrame);  
        x =  Ellipses{CurrentFrame}(CurrentEllipse,1)+1;
        y =  Ellipses{CurrentFrame}(CurrentEllipse,2)+1;
        if isempty(x)
            continue
        end 
        distances = ((px-x).^2 + (py-y).^2).^.5;
        candidate_indices = NucleusDistMat > distances; 
        %Record Fluorescence        
        NucleusIDMat(candidate_indices) = s;
        NucleusDistMat(candidate_indices) = distances(candidate_indices);
    end
    %Loop Through Particles to see which remain on and how long these have
    %been active
    ParticlesToShow = [];
    for i = 1:length(CompiledParticles)        
        cp_frames = CompiledParticles(i).Frame;        
        all_frames = min(cp_frames):max(cp_frames);    
        %the patch turns on when the spot is first spotted and then stays
        %on for the rest of the movie. 
        extant_frames = all_frames((all_frames <= CurrentFrame)&...
            (all_frames>first_frame)); 
        frame_filter = ismember(cp_frames,extant_frames);
        
        if any(frame_filter)               
            nucleus_filter = NucleusIDMat==CompiledParticles(i).Nucleus;            
            
            %this calculation could use some clarification
            if meanRate
                Fluo = nanmean(CompiledParticles(i).Fluo(frame_filter));                
            elseif instantRate
                if any(cp_frames==CurrentFrame)
                    Fluo = CompiledParticles(i).Fluo(cp_frames==CurrentFrame);
                else
                    Fluo = 0;
                end
            end
            Fluo = min(MaxmRNAFluo,Fluo) / MaxmRNAFluo;
            if fractionOn
                Fluo = .85;
            end
            
            for k = 1:3
                slice = NucleusStateMat(:,:,k);                        
                slice(nucleus_filter) = fill_color(k)*Fluo;
                NucleusStateMat(:,:,k) = slice;
            end                    
            ParticlesToShow = [ParticlesToShow CompiledParticles(i).Nucleus];            
        end
    end
    %Now Draw Nucleus Borders for active nuclei
    NucleusBorderMat = zeros(size(NucleusIDMat));
    window = 1; %radius of convolution window
    for i = ParticlesToShow
        %get coordinates of nucleus patch
        if sum(sum(NucleusIDMat==i)) > 0
            x_vec = reshape(px(NucleusIDMat==i),[],1);
            y_vec = reshape(py(NucleusIDMat==i),[],1);
            for j = 1:length(x_vec)
                neighborhoodScore = sum(sum(NucleusIDMat(max(1,y_vec(j)-window):min(y_vec(j)+window,yDim),...
                         max(1,x_vec(j)-window):min(x_vec(j) + window,xDim))));
                if neighborhoodScore~= i*(2*window+1)^2
                    NucleusBorderMat(y_vec(j),x_vec(j)) = 1;
                end
            end
        end        
    end
    %Prevent overlap between fluroescence mask and borders
    for k = 1:3
        slice = NucleusStateMat(:,:,k);
        slice(NucleusBorderMat>0) = 0;
        NucleusStateMat(:,:,k) = slice;
    end

    %Make a maximum projection of the mRNA channel
    D=dir([PreProcPath,filesep,Prefix,'_',num2str(CurrentFrame,'%03d'),'_z*.tif']);
    %Do not load the first and last frame as they are black
    ImageTemp=zeros(yDim, xDim, zSlices) ;
    for i=2:(length(D)-1)
        ImageTemp(:,:,i-1)=imread([PreProcPath,D(i).name]);
    end
    mRNAImageRaw=max(ImageTemp,[],3);    
    %Load the corresponding histone image
    HistoneImage=imread([PreProcPath,filesep,...
        Prefix,'-His_',num2str(CurrentFrame,'%03d'),'.tif']);        

    %Overlay all channels
    MCPshading = NucleusStateMat(:,:,2);  %this is the patch. should be dynamically scaled. 
    scale = max(MaxmRNAFluo) / 3;
    MCPChannel =  mRNAImageRaw/scale/2 + MCPshading; %this is the spot on top of the patch. 150 is the number that made the spot 
    %visible over the patch. 
    MCPChannel(MCPChannel>1) = 1;
         
    
    HistoneChannel=  mat2gray(HistoneImage)/2 + NucleusStateMat(:,:,1);%
    HistoneChannel(HistoneChannel>1) = 1;
    StateChannel = NucleusStateMat(:,:,3);    
    
    ImOverlay=cat(3,HistoneChannel,MCPChannel,StateChannel);
    imshow(fliplr(ImOverlay),'DisplayRange',[], 'Parent', OverlayAxes);
    ylim([16,yDim-16])
    if ceil(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))) == 40        
        saveas(gcf,[WritePath save_prefix '_patch_overlay_40min.tif']); 
    end
           
    text(OverlayAxes, 20,50,[num2str(round(ElapsedTime(CurrentFrame)-ElapsedTime(FrameRange(1))),'%02d'),...
        ' min'],'Color','k','FontSize',10,'BackgroundColor',[1,1,1,.5])	    %not sure about these numbers
    drawnow
    if saveFigs
        saveas(OverlayFig,[WritePath '\nc',num2str(nc),...
            '-',num2str(CurrentFrame,'%02d'),'.tif']);   
    end
    
end
