function paramsAn=masterfit3(mainfold,varargin)
% mainfold should have nd2 files or tiff stacks. if there are phasemask
% tiff files, they should correspond 1:1 with the nd2 files. the outputs
% written to disk are a 2Dtracks.fig, fits.fig, _analysis.mat, and a binary
% version of the nd2 (or tiff stack) movie.
%
% inside the _analysis.mat file there are the good fits, the tracked
% particles, the guessed particle positions, and a human-readable output
%
% there are various usage settings below.

%#ok<*PFBNS>

%% parameters

% CAMERA PARAMETERS
cParams.pixelSize = 49;
cParams.frameRate = 1/.04;

% ANALYSIS PARAMETERS
paramsAn.checkVals = 0;
paramsAn.phase = 1;
paramsAn.oldPhaseMasks = 1;
paramsAn.fitting = 0;
paramsAn.tracking = 1;

paramsAn.makeXls = 0;
paramsAn.bgsub = 0;
paramsAn.bgsubWidth = 50;
paramsAn.offset = 1000;
paramsAn.viewFits = inf;
paramsAn.minSep = 5;
paramsAn.widthLB = 1;
paramsAn.widthUB = 15;
paramsAn.loadOldParams = 0;
% fieldsAn = fieldnames(paramsAn);

% PHASEMASK PARAMETERS.
paramsPhase.nDilation = 1;
paramsPhase.lowThresh = 0;
paramsPhase.highThresh = 1;
paramsPhase.autofill = 1;
paramsPhase.minArea = 100;
paramsPhase.maxArea = 1e4;
fieldsPhase = fieldnames(paramsPhase);

% PEAK GUESSING PARAMETERS,
paramsPeaks.spotSizeLB = 1;
paramsPeaks.spotSizeUB = 10;
paramsPeaks.intThresh = 400;
paramsPeaks.hMax = 500;
paramsPeaks.lZero = 10;
fieldsPeaks = fieldnames(paramsPeaks);

% TRACKING PARAMETERS
paramsTr.minMerit = .08;
paramsTr.intTime = 1/cParams.frameRate;
paramsTr.gamma = 3;
paramsTr.minTrLength = 3;
paramsTr.maxStepSize = 10;
paramsTr.swh = 1;
paramsTr.delay = 0;
fieldsTr = fieldnames(paramsTr);

%% find all the movie files
display('Select the movies.')
[datalist,dataloc,~]=uigetfile([mainfold filesep '*.nd2;*.tif*;*.bin'],...
    'Matlab Files','multiselect','on');

if ~dataloc
    display('no data selected')
    return
end

if ~iscell(datalist); datalist={datalist}; end
for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
[dlocs,dnames,dexts]=cellfun(@fileparts,datalist,'uniformoutput',false);

%% WRITE BINARY FILES
if paramsAn.bgsub
    bFiles=dir([dataloc,'*_bgsub.bin']);
    for ii=1:numel(dlocs);
        if ~any(ismember({bFiles.name},[dnames{ii},'_bgsub.bin']))
            
            fid=fopen([fullfile(dlocs{ii},dnames{ii}) '_bgsub.bin'],'W');
            [~,nframes,sz]=binGetFrames2([fullfile(dlocs{ii},dnames{ii}),'.bin'],1);
            
            for jj=1:floor(nframes/subwidth)
                waitbar(jj/floor(nframes/subwidth),h1)
                if jj~=floor(nframes/subwidth)
                    v=binGetFrames2([fullfile(dlocs{ii},dnames{ii}) '.bin'],...
                        1+subwidth*(jj-1):subwidth*jj);
                else
                    v=binGetFrames2([fullfile(dlocs{ii},dnames{ii}) '.bin'],...
                        1+subwidth*(jj-1):nframes);
                end
                
                fwrite(fid,uint16(bsxfun(@plus,double(v),-mean(v,3))+offset),'uint16');
            end
            fwrite(fid,sz,'double');
            fwrite(fid,nframes,'double');
            
            fclose(fid);
            
            % rename the data file names
            dnames{ii}=[dnames{ii},'_bgsub'];
        end
        
        % rename the data names variable
        dnames{ii}=[dnames{ii},'_bgsub'];
    end
elseif paramsAn.bgsub == 0
    bFiles=dir([dataloc,'*.bin']);
    for ii=1:numel(dnames);
        if ~any(ismember({bFiles.name},[dnames{ii},'.bin']))
            writebin(fullfile(dlocs{ii},[dnames{ii},dexts{ii}]));
        end
    end
end

%% initialize analysis files
if paramsAn.phase
    display('Select the phase images.')
    [phaselist,phaselistloc,~]=uigetfile([mainfold filesep...
        '*.nd2;*.tif*;*.bin'],'Matlab Files','multiselect','on');
    
    if ~phaselistloc
        display('no data selected')
        return
    end
    
    if ~iscell(phaselist); phaselist={phaselist}; end
    for ii=1:numel(phaselist); phaselist{ii}=[phaselistloc phaselist{ii}]; end
    [plocs,pnames,pexts]=cellfun(@fileparts,phaselist,'uniformoutput',false);
else
    for ii=1:numel(dnames)
        [~,~,sz]=binGetFrames2([fullfile(dlocs{ii},dnames{ii}),'.bin'],1);
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        phaseMask=ones(sz);
        m.phaseMask=phaseMask;
    end
end

%% compute phase mask
if paramsAn.phase&&~paramsAn.oldPhaseMasks
    for ii=1:numel(dnames)
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        if strcmp(pexts{ii},'.nd2')
            vidid=bfGetReader(phaselist{ii});
            img=bfGetPlane(vidid,ii);
        else
            img=imread(fullfile(plocs{ii},[pnames{ii},pexts{ii}]));
        end
        
        going = 1;
        while going
            oldParams = [paramsPhase.nDilation,paramsPhase.lowThresh,paramsPhase.highThresh,...
                paramsPhase.autofill,paramsPhase.minArea,paramsPhase.maxArea];
            
            phaseMask=valley2(img,paramsPhase);
            
            if paramsAn.checkVals
                subplot(121)
                imshow(phaseMask,[])
                subplot(122)
                imshow(img,[])
                
                % prompt user for parameter changes
                dlgPrompt=fieldsPhase;
                dlgTitle='Phasemask Parameters';
                numDlgLines=1;
                dValues = cellfun(@num2str,num2cell(oldParams),'uniformoutput',0);
                opts.WindowStyle='normal';
                
                newParams=cellfun(@str2double,...
                    inputdlg(dlgPrompt,dlgTitle,numDlgLines,dValues,opts))';
                
                if any(newParams~=oldParams)
                    paramsPhase.nDilation = newParams(1);
                    paramsPhase.lowThresh = newParams(2);
                    paramsPhase.highThresh = newParams(3);
                    paramsPhase.autofill = newParams(4);
                    paramsPhase.minArea = newParams(5);
                    paramsPhase.maxArea = newParams(6);
                else
                    going = 0;      % advance to next phase image if no parameters changed
                end
            end
        end
        
        [~,goodcells]=select_cells(phaseMask,img);
        phaseMask(logical(-1*(ismember(phaseMask,goodcells)-1)))=0;
        
        m.phaseMask=phaseMask;
        m.paramsPhase = paramsPhase;
    end
end

%% analyze movie files
for currMovie=1:numel(dnames)
    [~,nFrames,~]=binGetFrames2([fullfile(...
        dlocs{currMovie},dnames{currMovie}),'.bin'],1);
    
    m=matfile([fullfile(dlocs{currMovie},dnames{currMovie}),...
        '_analysis.mat'],'Writable',true);
    phaseMask=m.phaseMask;
    phaseImg=imread(fullfile(plocs{currMovie},[pnames{currMovie},pexts{currMovie}]));
    
    goodFits = zeros(0,23);
    if  ~paramsAn.fitting
        try
            goodFits = m.goodFits;
            guesses = m.guesses;
        catch
            warning(['missing data in analysis file number ' num2str(currMovie)])
            continue
        end
    else
        %% fitting
        fprintf(['Fitting this movie: ',dnames{currMovie},'\n'])
        
        p = cell(1,nFrames);
        ers = cell(1,nFrames);
        guesses = cell(1,nFrames);
        if paramsAn.checkVals
            % parameter checking version
            currFrame = 1;
            while currFrame<=nFrames
                thisframe=binGetFrames2(...
                    [fullfile(dlocs{currMovie},dnames{currMovie}) '.bin'],currFrame);
                
                [p{currFrame},ers{currFrame},guesses{currFrame}]=gaussFit(thisframe,phaseMask,...
                    'spotSizeLB',paramsPeaks.spotSizeLB,...
                    'spotSizeUB',paramsPeaks.spotSizeUB,...
                    'intThresh',paramsPeaks.intThresh,...
                    'hMax',paramsPeaks.hMax,...
                    'lZero',paramsPeaks.lZero,...
                    'showGuessing',paramsAn.checkVals,...
                    'frameNumber',currFrame);
                
                % prompt user for parameter changes (and advance to next frame)
                oldParams = [paramsPeaks.spotSizeLB,paramsPeaks.spotSizeUB,...
                    paramsPeaks.intThresh,paramsPeaks.hMax,paramsPeaks.lZero,currFrame+1];
                dlgPrompt=[fieldsPeaks;'next frame number'];
                dlgTitle='Peak Guessing. Close to advance to next movie.';
                numDlgLines=1;
                dValues=cellfun(@num2str,num2cell(oldParams),'uniformoutput',0);
                opts.WindowStyle='normal';
                
                newParams = cellfun(@str2double,...
                    inputdlg(dlgPrompt,dlgTitle,numDlgLines,dValues,opts))';
                
                if ~isempty(newParams)&&any(newParams ~= oldParams)
                    paramsPeaks.spotSizeLB = newParams(1);
                    paramsPeaks.spotSizeUB = newParams(2);
                    paramsPeaks.intThresh = newParams(3);
                    paramsPeaks.hMax = newParams(4);
                    paramsPeaks.lZero = newParams(5);
                    newParams(6) = newParams(6) - 1;    % don't advance to next frame
                end
                
                if isempty(newParams)
                    currFrame = nFrames+1;              % advance to next movie
                else
                    currFrame = newParams(6);           % advance to user input frame.
                end
            end
            
        elseif ~paramsAn.checkVals
            % parfor loop version
            parfor currFrame = 1:nFrames
                thisframe=binGetFrames2(...
                    [fullfile(dlocs{currMovie},dnames{currMovie}) '.bin'],currFrame);
                
                [p{currFrame},ers{currFrame},guesses{currFrame}]=gaussFit(thisframe,phaseMask,...
                    'spotSizeLB',paramsPeaks.spotSizeLB,...
                    'spotSizeUB',paramsPeaks.spotSizeUB,...
                    'intThresh',paramsPeaks.intThresh,...
                    'hMax',paramsPeaks.hMax,...
                    'lZero',paramsPeaks.lZero,...
                    'showGuessing',paramsAn.checkVals,...
                    'frameNumber',currFrame);
            end
            
            % results selection and reorganization to legacy format
            orgFun = @(in_p,in_ers,nf,fn,pm) cat(2,...
                fn*ones(nf,1),...
                (1:size(in_p,1))',...
                in_p(:,1),in_ers(:,1),...
                in_p(:,2),in_ers(:,2),...
                in_p(:,3),in_ers(:,3),...
                in_p(:,4),in_ers(:,4),...
                in_p(:,5),in_ers(:,5),...
                nan(nf,1),nan(nf,1),nan(nf,1),...
                nan(nf,1),nan(nf,1),...
                nan(nf,1),nan(nf,1),nan(nf,1),...
                pm(sub2ind(size(pm),ceil(in_p(:,2)),ceil(in_p(:,1)))),...
                nan(nf,1),nan(nf,1));
            
            allFits = cellfun(@(x,y,z,w)orgFun(x,y,z,w,phaseMask),...
                p,ers,num2cell(cellfun(@(x)size(x,1),p)),num2cell(1:nFrames),...
                'uniformoutput',0);
            allFits = cat(1,allFits{:});
            
            whichGood = ...
                allFits(:,7) > paramsAn.widthLB & ...
                allFits(:,7) < paramsAn.widthUB & ...
                allFits(:,21) > 0;
            
            goodFits = allFits(whichGood,:);
            
            % WRITE ANALYSIS FILE
            m.goodFits = goodFits;
            m.allFits = allFits;
            m.guesses = guesses;
            
            imshow(phaseImg,[]); hold all
            scatter(goodFits(:,3),goodFits(:,5),'.'); hold off
            saveas(gcf,[fullfile(dlocs{currMovie},dnames{currMovie}) '_goodFits.fig'])
        end
        
        % record the changes to the peak finding parameters
        m.paramsPeaks = paramsPeaks;
    end
    
    %% Tracking
    if ~paramsAn.tracking
        try
            trfile = m.trackfile;
        catch
            trfile = zeros(0,16);
            warning(['no tracking data in analysis file number ' num2str(currMovie)])
            continue
        end
    else
        if paramsAn.tracking
            % if goodfits is empty, don't bother tracking.
            if isempty(goodFits)
                going = 0;
            else
                going = 1;
            end
            
            while going
                oldParams = [paramsTr.minMerit, paramsTr.intTime, paramsTr.gamma, ...
                    paramsTr.minTrLength, paramsTr.maxStepSize, paramsTr.swh, ...
                    paramsTr.delay];
                alpha=-log(oldParams(1))/oldParams(5);
                trfile=Track_3D2(goodFits,[],[],oldParams(1),...
                    alpha,oldParams(3),oldParams(4),oldParams(6),...
                    cParams.pixelSize,oldParams(7),oldParams(2));
                
                if isempty(trfile)
                    fprintf(['No available tracks for ''',dnames{currMovie}...
                        '''. Check tracking parameters.\n']);
                else
                    hastrack=unique(trfile(:,1))';
                    cm=jet(numel(hastrack));
                    
                    imshow(phaseImg,[]); hold all
                    
                    for trnum = hastrack
                        plot(trfile(trfile(:,1)==trnum,4),...
                            trfile(trfile(:,1)==trnum,5),...
                            'Color',cm(hastrack==trnum,:),...
                            'linewidth',1);
                    end
                    title([num2str(numel(hastrack)) ' tracks'])
                    hold off
                    
                    saveas(gcf,[fullfile(dlocs{currMovie},dnames{currMovie}),'_tracks.fig']);
                end
                
                if ~paramsAn.checkVals
                    going = 0;      % progress to next movie if parameters are assumed correct
                else
                    % prompt user for parameter changes
                    dlgPrompt=fieldsTr;
                    dlgTitle='Tracking';
                    numDlgLines=[1,length(dlgTitle)+10];
                    dValues={num2str(oldParams(1)),num2str(oldParams(2)),...
                        num2str(oldParams(3)),num2str(oldParams(4)),...
                        num2str(oldParams(5)),num2str(oldParams(6)),...
                        num2str(oldParams(7))};
                    opts.WindowStyle='normal';
                    
                    newParams = cellfun(@str2double,inputdlg(dlgPrompt,...
                        dlgTitle,numDlgLines,dValues,opts))';
                    if isempty(newParams)
                        going = 0;
                    end
                    
                    if ~isempty(newParams)&&any(newParams~=oldParams)
                        paramsTr.minMerit = newParams(1);
                        paramsTr.intTime = newParams(2);
                        paramsTr.gamma = newParams(3);
                        paramsTr.minTrLength = newParams(4);
                        paramsTr.maxStepSize = newParams(5);
                        paramsTr.swh = newParams(6);
                        paramsTr.delay = newParams(7);
                    else
                        going = 0;  % progress to next movie if parameters are left unchanged
                    end
                end
                
                m.paramsTr=paramsTr;
            end
            
            % record tracking results
            m.trackfile=trfile;
        end
    end
    
    %% misc. outputs
    % Output ViewFit Files for the Current Movie (If Selected)
    if any(ismember(paramsAn.viewFits,currMovie))||isinf(paramsAn.viewFits)
        % Output ViewFits frames if the fit file is not empty
        display(['Generating ViewFits frames for ''',...
            fullfile(dlocs{currMovie},dnames{currMovie}),'''...'])
        
        if paramsAn.phase
            bounds = [find(any(phaseMask,2),1,'first'),find(any(phaseMask,2),1,'last'),...
                find(any(phaseMask,1),1,'first'),find(any(phaseMask,1),1,'last')];
        else
            bounds = [1,size(phaseMask,1),1,size(phaseMask,2)];
        end
        
        Viewfits3([fullfile(dlocs{currMovie},dnames{currMovie}),'.bin'],...
            goodFits,guesses,trfile,cParams.frameRate,bounds);
    end
    
    if paramsAn.makeXls
        tstr=[{'Frame'} {'Molecule'} {'Amplitude'} {'+/-'} {'Offset'}...
            {'+/-'} {'X-Width'} {'+/-'} {'X Center'} {'+/-'}  {'Y Center'}...
            {'+/-'} {'Good Fit?'} {'Integral'} {'Small box width'}...
            {'Y-Width'} {'+/-'} {'Wx/Wy'} {'Z Center'} {'+/-'} {'Cell No'}...
            {'sROI'} {'tROI'}];% {'distance from center'}];
        xlswrite([dnames{currMovie}],cat(1,tstr,num2cell(goodFits)))
    end
end
end

function Viewfits3(vidloc,fitfile1,guesses,trfile,framerate,bounds)
% This code takes in raw single-molecule tif frames and fits files and put
% a square at each fit. Useful for checking to see if fitting parameters
% are right.

% INPUTS:

% OUTPUTS:

%-------------------------------------------------------------------------%
    function cFun(x,y)
        close(x)
        close(y)
    end

halfBoxWidth = 7;
intMax=1;

[~,nframes,~]=binGetFrames2(vidloc,1);
sz = [bounds(2)-bounds(1)+1,bounds(4)-bounds(3)+1];

wobj=VideoWriter([vidloc(1:end-4),'_Viewfits'],'Uncompressed AVI');
wobj.FrameRate=framerate;
% wobj.Quality=100;

open(wobj)
h=waitbar(0,['viewfits for movie ' vidloc(1:end-4)]);
c=onCleanup(@()cFun(h,wobj));
for ii=1:nframes
    if rem(ii,10)==0
        waitbar(ii/nframes,h)
    end
    
    img=double(binGetFrames2(vidloc,ii));
    img=img(bounds(1):bounds(2),bounds(3):bounds(4));
    
    % autoscale
    img=img-min(img(:));
    img=img/max(img(:));
    
    % Loop for all guess peaks in the current frame
    for jj = 1:size(guesses{ii},1)
        pix_x=guesses{ii}(jj,2)-bounds(3);
        pix_y=guesses{ii}(jj,1)-bounds(1);
        if pix_x<sz(2)&&pix_y<sz(1)&&pix_x>0&&pix_y>0
            img=makebox(img,intMax/2,pix_x,pix_y,halfBoxWidth,sz);
        end
    end
    
    if any(size(img)~=sz)
        keyboard
    end
    
    % Indexing which rows in the fits file belong to the current frame
    r=find(fitfile1(:,1)==ii);
    
    % Loop for all good fits in the current frame
    for jj=1:numel(r)
        pix_x=round(fitfile1(r(jj),5))-bounds(3);
        pix_y=round(fitfile1(r(jj),3))-bounds(1);
        
        img=makebox(img,intMax,pix_x,pix_y,halfBoxWidth,sz);
    end
    
    if any(size(img)~=sz)
        keyboard
    end
    
    if numel(trfile)>0
        tr=trfile(trfile(:,2)==ii,:);
        
        cm=jet(5);
        img=img(:,:,[1,1,1]);
        for jj=1:size(tr,1)
            pix_x=round(tr(jj,5))-bounds(3);
            pix_y=round(tr(jj,4))-bounds(1);
            img=makebox(img,cm(mod(tr(jj,1),size(cm,1)-1)+1,:),...
                pix_x,pix_y,halfBoxWidth,sz);
        end
    end
    
    if any(size(img(:,:,1))~=sz)
        keyboard
    end
    
    %     imwrite(currentframe2, [vidloc(1:end-4),'_Viewfits.tif'], 'writemode', 'append');
    writeVideo(wobj,img)
end
close(wobj);
end

function currentframe=makebox(currentframe,curr_fr_maxint,...
    pix_x,pix_y,half_symbol_size,sz)

for ii=1:size(currentframe,3)
    b=pix_y-half_symbol_size;
    t=pix_y+half_symbol_size;
    l=pix_x-half_symbol_size;
    r=pix_x+half_symbol_size;
    
    b(b<1)=1;
    l(l<1)=1;
    t(t>sz(1))=sz(1);
    r(r>sz(2))=sz(2);
    
    currentframe(b:t,[l,r],ii)=curr_fr_maxint(ii);
    currentframe([b,t],l:r,ii)=curr_fr_maxint(ii);
end
end

function [cell_xy,good_cell]=select_cells(PhaseMask,img)

% Let user click on the phase mask of cell images to decide which cell to
% analyze subsequently

% INPUTS:

% PhaseMask: Numbered phase mask of cell images returned by the 'valley.m'
% code or any other cell segmentation code

% OUTPUTS:

% cell_xy: The x- and y-coordinates clicked on by the users to select cells

% good_cell: The indices of chosen cells

cell_fig_h=figure;
subplot(121);
imshow(img,[],'Border','Tight');
subplot(122);
imshow(PhaseMask~=0,[],'Border','Tight');
title(sprintf(['Left-click on cells to be analyzed.\n Press ''Enter'' to'...
    ' proceed.\n Or click return without clicking\n on any cell to analyze ALL of them.']))

% Now move the axes slightly so that the top of the title is visible
set(cell_fig_h,'Units','normalized')
P=get(cell_fig_h,'Position');
set(cell_fig_h,'Position',P+[0 -0.03 0 0])

conv_hull=regionprops(PhaseMask,'Convexhull');
for i=1:length(conv_hull)
    hold all,
    plot(conv_hull(i,1).ConvexHull(:,1),...
        conv_hull(i,1).ConvexHull(:,2),'c-','linewidth',2)
end

key_input=0;
% Let the user pick cells with good shapes
cell_x=[]; cell_y=[];
while key_input~=121
    [curr_cell_x,curr_cell_y,key_input]=ginput(1);
    hold all,
    plot(curr_cell_x,curr_cell_y,'m*','markersize',12);
    cell_x=vertcat(cell_x,curr_cell_x);  %#ok<AGROW>
    cell_y=vertcat(cell_y,curr_cell_y); %#ok<AGROW>
end

if ~isempty(cell_x) % If the user did click on something
    cell_xy=[cell_x,cell_y];
    
    good_cell=PhaseMask(sub2ind(size(PhaseMask),round(cell_y),round(cell_x)));
    good_cell=good_cell(good_cell~=0);
    if ~isempty(good_cell) % If the user did click on at least 1 cell
        good_cell=unique(good_cell);
    else % If the user only clicked on background, use all cells.
        good_cell=1:max(PhaseMask(:));
    end
else % If the user hits "Enter" directly, use all cells.
    cell_xy=[];
    good_cell=1:max(PhaseMask(:));
end
close(cell_fig_h);
end