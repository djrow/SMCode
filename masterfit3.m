function paramsAn=masterfit3(mainfold)
% mainfold should have nd2 files or tiff stacks. if there are phasemask
% tiff files, they should correspond 1:1 with the nd2 files. the outputs
% written to disk are a 2Dtracks.fig, fits.fig, _analysis.mat, and a binary
% version of the nd2 (or tiff stack) movie.
%
% inside the _analysis.mat file there are the good fits, the tracked
% particles, the guessed particle positions, and a human-readable output
%
% there are various usage settings below.

%% parameters
% CAMERA PARAMETERS
cParams.pixelSize = 49;
cParams.frameRate = 1/.01;

% ANALYSIS PARAMETERS
paramsAn.checkVals = 0;
paramsAn.phaseImages = 0;
paramsAn.okWithPhase = 1;
paramsAn.fittingBool = 1;
paramsAn.trackingBool = 1;
paramsAn.makeXls = 0;
paramsAn.viewFits = 'all'; % 'all';

paramsAn.bgsub = 0;
paramsAn.bgsubWidth = 50;
paramsAn.offset = 1000;

paramsAn.widthLB = 1;
paramsAn.widthUB = 15;

% PHASEMASK PARAMETERS.
paramsPhase.nDilation = 1;
paramsPhase.lowThresh = 0;
paramsPhase.highThresh = 1;
paramsPhase.autofill = 1;
paramsPhase.minArea = 100;
paramsPhase.maxArea = 1e4;
fieldsPhase = fieldnames(paramsPhase);

% PEAK GUESSING PARAMETERS,
paramsPeaks.spotSizeLB = 1.2;
paramsPeaks.spotSizeUB = 8;
paramsPeaks.intThresh = 500;
paramsPeaks.hMax = 700;
paramsPeaks.lZero = 10;
fieldsPeaks = fieldnames(paramsPeaks);

% TRACKING PARAMETERS
paramsTr.minMerit = .1;
paramsTr.intTime = 1/cParams.frameRate;
paramsTr.gamma = .25;
paramsTr.minTrLength = 40;
paramsTr.maxStepSize = 20;
paramsTr.swh = 1;
paramsTr.delay = 0;
fieldsTr = fieldnames(paramsTr);

% reorganization to legacy format
orgFun = @(in_p,in_ers,nf,fn,pIDs,gs) cat(2,...
    fn*ones(nf,1),(1:nf)',...       % frame number, fit number
    in_p(:,1),in_ers(:,1),...       % center 1, center 1 error
    in_p(:,2),in_ers(:,2),...       % center 2, center 2 error
    in_p(:,3),in_ers(:,3),...       % s.d., s.d. error
    in_p(:,4),in_ers(:,4),...       % amplitude, amplitude error
    in_p(:,5),in_ers(:,5),...       % offset, offset error
    gs(:,1),gs(:,2),...             % guessesX, guessesY
    nan(nf,1),nan(nf,1),nan(nf,1),...
    nan(nf,1),nan(nf,1),nan(nf,1),... % 6 empty columns
    pIDs, ...                       % cell number
    nan(nf,1),nan(nf,1));           % empty columns

% pay no attention to the man behind the curtain
%#ok<*PFBNS>
%#ok<*NASGU>
d = onCleanup(@()eval('fclose all;'));
opts.WindowStyle = 'normal';
numDlgLines = 1;

%% find all the movie files
[datalist,dataloc,~]=uigetfile([mainfold filesep '*.nd2;*.tif*;*.bin'],...
    'Select the movies','multiselect','on');

if ~dataloc
    display('no data selected')
    return
end

if ~iscell(datalist); datalist={datalist}; end
datalist = cellfun(@(x)[dataloc, x],datalist,'uniformoutput',false);
[dlocs,dnames,dexts]=cellfun(@fileparts,datalist,'uniformoutput',false);

%% initialize analysis *.mat files
if paramsAn.phaseImages
    [phaselist,phaselistloc,~]=uigetfile([mainfold filesep...
        '*.nd2;*.tif*;*.bin'],'Select the phase images','multiselect','on');
    
    if ~phaselistloc
        display('no data selected')
        return
    end
    
    if ~iscell(phaselist); phaselist={phaselist}; end
    phaselist = cellfun(@(x)[phaselistloc, x],phaselist,'uniformoutput',false);
    [plocs,pnames,pexts]=cellfun(@fileparts,phaselist,'uniformoutput',false);
    
    for ii=1:numel(pnames)
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        if strcmp(pexts{ii},'.nd2')
            vidid = bfGetReader(phaselist{ii});
            phaseImg = bfGetPlane(vidid,ii);
        else
            phaseImg=imread(fullfile(plocs{ii},[pnames{ii},pexts{ii}]));
        end
        
        m.phaseImg = phaseImg;
    end
% else
%     for ii=1:numel(dnames)
% %         [~,~,sz]=binGetFrames2([fullfile(dlocs{ii},dnames{ii}),'.bin'],1);
%         m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
%         
%         m.phaseMask = ones(128);
%         m.phaseImg = ones(128);
%     end
end

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
        end
        
        % rename the data names variable
        dnames{ii}=[dnames{ii},'_bgsub'];
    end
else
    bFiles=dir([dataloc,'*.bin']);
    for ii=1:numel(dnames);
        if ~any(ismember({bFiles.name},[dnames{ii},'.bin']))
            writebin(fullfile(dlocs{ii},[dnames{ii},dexts{ii}]));
        end
    end
end

%% compute phase mask
if paramsAn.phaseImages&&~paramsAn.okWithPhase
    for ii=1:numel(dnames)
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        phaseImg=m.phaseImg;
        
        going = 1;
        while going
            oldParams = [paramsPhase.nDilation,paramsPhase.lowThresh,paramsPhase.highThresh,...
                paramsPhase.autofill,paramsPhase.minArea,paramsPhase.maxArea];
            
            phaseMask = valley2(phaseImg,paramsPhase);
            
            subplot(2,4,1)
            imshow(phaseMask,[])
            subplot(2,4,5)
            imshow(phaseImg,[])
            
            % prompt user for parameter changes
            dValues = cellfun(@num2str,num2cell(oldParams),'uniformoutput',0);
            newParams=cellfun(@str2double,...
                inputdlg(fieldsPhase,'Phasemask Parameters',...
                numDlgLines,dValues,opts))';
            
            if isempty(newParams)
                % advance to next phase image if box is closed
                going = 0;
            elseif any(newParams~=oldParams)
                % change parameters if parameters changed
                paramsPhase.nDilation = newParams(1);
                paramsPhase.lowThresh = newParams(2);
                paramsPhase.highThresh = newParams(3);
                paramsPhase.autofill = newParams(4);
                paramsPhase.minArea = newParams(5);
                paramsPhase.maxArea = newParams(6);
            else
                % advance to next phase image if no parameters changed
                going = 0;
            end
        end
        
        phaseMask = selectCells(phaseMask,phaseImg);
        
        m.phaseMask = phaseMask;
        m.paramsPhase = paramsPhase;
    end
end

%% analyze movie files
for movieNum = 1:numel(dnames)
    % collect information from disk
    movieName = fullfile(dlocs{movieNum},dnames{movieNum});
    [~,nFrames,vSize] = binGetFrames2([movieName,'.bin'],1);
    m = matfile([movieName, '_analysis.mat'],'Writable',true);
    
    
    try
        phaseMask = m.phaseMask;
        phaseImg = m.phaseImg;
    catch
        phaseMask = ones(vSize);
        phaseImg = ones(vSize);
        m.phaseMask = phaseMask;
        m.phaseImg = phaseMask;
    end
    
    if paramsAn.fittingBool
        allFits = zeros(0,23);
        goodFits = zeros(0,23);
        %% fitting
        fprintf(['Fitting this movie: ',dnames{movieNum},'\n'])
        
        % a small multicellular organism is initialized
        p = cell(1,nFrames); p(:) = {nan(1,5)}; guesses = p;
        outPut = p; allFitsCell = p; nFits = p;
        if paramsAn.checkVals
            currFrame = 1;
            
            dlgPrompt = [fieldsPeaks; 'next frame number'];
            while currFrame <= nFrames
                thisframe=double(binGetFrames2([movieName '.bin'],currFrame));
                
                [~] = gaussFit(thisframe,...
                    'spotSizeLB',paramsPeaks.spotSizeLB,...
                    'spotSizeUB',paramsPeaks.spotSizeUB,...
                    'intThresh',paramsPeaks.intThresh,...
                    'hMax',paramsPeaks.hMax,...
                    'lZero',paramsPeaks.lZero,...
                    'checkVals',paramsAn.checkVals,...
                    'searchBool',1,...
                    'frameNumber',currFrame);
                
                % prompt user for parameter changes
                oldParams = [paramsPeaks.spotSizeLB,paramsPeaks.spotSizeUB,...
                    paramsPeaks.intThresh,paramsPeaks.hMax,paramsPeaks.lZero,currFrame+1];
                dValues = cellfun(@num2str,num2cell(oldParams),'uniformoutput',0);
                newParams = cellfun(@str2double,...
                    inputdlg(dlgPrompt,'Peak Guessing',numDlgLines,dValues,opts))';
                
                if ~isempty(newParams)&&any(newParams ~= oldParams)
                    % don't advance to next frame
                    paramsPeaks.spotSizeLB = newParams(1);
                    paramsPeaks.spotSizeUB = newParams(2);
                    paramsPeaks.intThresh = newParams(3);
                    paramsPeaks.hMax = newParams(4);
                    paramsPeaks.lZero = newParams(5);
                    newParams(6) = newParams(6) - 1;
                end
                
                if isempty(newParams)
                    % advance to next movie if box is closed
                    currFrame = nFrames+1;
                else
                    % advance to user input frame.
                    currFrame = newParams(6);
                end
            end
        elseif ~paramsAn.checkVals
            % parfor loop version
            parfor currFrame = 1:nFrames
                thisframe = double(binGetFrames2([movieName '.bin'],currFrame));
                
                [p{currFrame},guesses{currFrame}] = ...
                    gaussFit(thisframe,...
                    'spotSizeLB',paramsPeaks.spotSizeLB,...
                    'spotSizeUB',paramsPeaks.spotSizeUB,...
                    'intThresh',paramsPeaks.intThresh,...
                    'hMax',paramsPeaks.hMax,...
                    'lZero',paramsPeaks.lZero,...
                    'checkVals',paramsAn.checkVals,...
                    'searchBool',1,...
                    'frameNumber',currFrame);
                
%                 [p{currFrame},guesses{currFrame}] = gaussFit(thisframe);
            end
            display('done fitting')
            
            % cell id assignments
            pmID = cellfun(@(x)findPmID(phaseMask,vSize,x),p, 'uniformoutput',false);
            
            nFits = cellfun(@(x)size(x,1),p, 'uniformoutput',false);
            ers = p;
            
            % reorganized outputs
            allFitsCell = cellfun(orgFun,...
                p,ers,nFits,num2cell(1:nFrames),pmID,guesses, 'uniformoutput',0);
            allFits = cat(1,allFitsCell{:});
            
            % results selection
            whichGood = ...
                allFits(:,7) > paramsAn.widthLB & ...   % fit is within the width tolerance
                allFits(:,7) < paramsAn.widthUB & ...   % fit is within the width tolerance
                allFits(:,21) > 0;                      % fit is in a cell
            goodFits = allFits(whichGood,:);
            
            % WRITE ANALYSIS FILE
            m.goodFits = goodFits;
            m.allFits = allFits;
            %             m.outPut = outPut;
        end
        
        m.paramsPeaks = paramsPeaks;
    elseif ~paramsAn.fittingBool
        try
            allFits = m.allFits;
            goodFits = m.goodFits;
        catch
            allFits = zeros(0,23);
            goodFits = zeros(0,23);
            warning(['missing data in analysis file number ' num2str(movieNum)])
            continue
        end
    end
    
    %% tracking
    if paramsAn.trackingBool
        if isempty(goodFits)
            going = 0;
            trackfile = zeros(0,16);
            warning('goodfits file is empty, but you wished to track')
        else
            going = 1;
        end
        
        while going
            oldParams = [paramsTr.minMerit, paramsTr.intTime, paramsTr.gamma, ...
                paramsTr.minTrLength, paramsTr.maxStepSize, paramsTr.swh, ...
                paramsTr.delay];
            alpha=-log(oldParams(1))/oldParams(5);
            trackfile=Track_3D2(goodFits,[],[],oldParams(1),...
                alpha,oldParams(3),oldParams(4),oldParams(6),...
                cParams.pixelSize,oldParams(7),oldParams(2));
            
            if isempty(trackfile)
                warning(['No tracks found in ''', movieName ...
                    ']''. Check tracking parameters'])
            else
                hastrack=unique(trackfile(:,1))';
                cm=jet(numel(hastrack));
                
                imshow(phaseImg,[]); hold all
                
                for trnum = hastrack
                    plot(trackfile(trackfile(:,1)==trnum,5),trackfile(trackfile(:,1)==trnum,4),...
                        'Color',cm(hastrack==trnum,:),'linewidth',1);
                end
                title([num2str(numel(hastrack)) ' tracks'])
                hold off
            end
            
            if paramsAn.checkVals
                % prompt user for parameter changes
                dValues = {num2str(oldParams(1)),num2str(oldParams(2)),...
                    num2str(oldParams(3)),num2str(oldParams(4)),...
                    num2str(oldParams(5)),num2str(oldParams(6)),...
                    num2str(oldParams(7))};
                newParams = cellfun(@str2double,inputdlg(fieldsTr,...
                    'Tracking',numDlgLines,dValues,opts))';
                if isempty(newParams)
                    going = 0;
                end
                
                if ~isempty(newParams) && any(newParams ~= oldParams)
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
            else
                % progress to next movie if parameters are assumed correct
                going = 0;
            end
            
            saveas(gcf,[fullfile(dlocs{movieNum},dnames{movieNum}),'_tracks.fig']);
            m.paramsTr = paramsTr;
        end
        
        % record tracking results
        m.trackfile = trackfile;
    elseif ~paramsAn.trackingBool
        try
            trackfile = m.trackfile;
        catch
            trackfile = zeros(0,16);
            warning(['no tracking data in analysis file for ' movieName(end-10:end)])
            continue
        end
    end
    
    %% misc. outputs
    % Output ViewFit Files for the Current Movie (If Selected)
    if any(ismember(paramsAn.viewFits,movieNum))||ischar(paramsAn.viewFits)
        % Output ViewFits frames if the fit file is not empty
        display(['Writing ViewFits movie for ''', movieName ,'''...'])
        
        if paramsAn.phaseImages
            bounds = [find(any(phaseMask,2),1,'first'),find(any(phaseMask,2),1,'last'),...
                find(any(phaseMask,1),1,'first'),find(any(phaseMask,1),1,'last')];
        else
            bounds = [1,size(phaseMask,1),1,size(phaseMask,2)];
        end
        
        Viewfits3([fullfile(dlocs{movieNum},dnames{movieNum}),'.bin'],...
            allFits,goodFits,trackfile,cParams.frameRate,bounds);
    end
    
    if paramsAn.makeXls
        tstr=[{'Frame'} {'Molecule'} {'Amplitude'} {'+/-'} {'Offset'}...
            {'+/-'} {'X-Width'} {'+/-'} {'X Center'} {'+/-'}  {'Y Center'}...
            {'+/-'} {'Good Fit?'} {'Integral'} {'Small box width'}...
            {'Y-Width'} {'+/-'} {'Wx/Wy'} {'Z Center'} {'+/-'} {'Cell No'}...
            {'sROI'} {'tROI'}];
        xlswrite([dnames{movieNum}],cat(1,tstr,num2cell(goodFits)))
    end
end
end

%% AUXILIARY FUNCTIONS

function Viewfits3(vidloc,allFits,goodFits,trackFile,framerate,bounds)
% This code takes in raw single-molecule tif frames and fits files and put
% a square at each fit. Useful for checking to see if fitting parameters
% are right.

% INPUTS:

% OUTPUTS:

colNums.guesses = [13,14];
colNums.fits = [3,5];
colNums.tracks = [4,5];

halfBoxWidth = 7;
intMax = 1;
cm = jet(5);

[~, nframes, ~] = binGetFrames2(vidloc, 1);
sz = [bounds(2)-bounds(1)+1, bounds(4)-bounds(3)+1];

wobj=VideoWriter([vidloc(1:end-4),'_Viewfits'],'Uncompressed AVI');
wobj.FrameRate=framerate;

open(wobj)
d = onCleanup(@()close(wobj));

h = waitbar(0,'writing viewfits');
c = onCleanup(@()close(h));

for ii = 1:nframes
    if rem(ii,10) == 0
        waitbar(ii/nframes,h,'writing viewfits')
    end
    
    % which rows in the fits files belong to the current frame
    r1 = find(allFits(:,1)==ii);        % guesses are stored in the allFits file
    r2 = find(goodFits(:,1)==ii);       % good fits
    r3 = find(trackFile(:,2)==ii);      % tracks, second column is frame number
    
    % guesses, fits, and tracks from current frame
    guessesX = allFits(r1,colNums.guesses(1)) - bounds(1);
    guessesY = allFits(r1,colNums.guesses(2)) - bounds(3);
    gfX = round(goodFits(r2,colNums.fits(1)) - bounds(1));
    gfY = round(goodFits(r2,colNums.fits(2)) - bounds(3));
    trX = round(trackFile(r3,colNums.tracks(1)) - bounds(1));
    trY = round(trackFile(r3,colNums.tracks(2)) - bounds(3));
    
    % load image frame from movie
    img = double(binGetFrames2(vidloc,ii));
    img = img(bounds(1):bounds(2),bounds(3):bounds(4));
    
    % autoscale
    img = img-min(img(:));
    img = img/max(img(:));
    
    % Loop for all guess peaks in the current frame
    for jj = 1:numel(r1)
        img = makeBox(img, 1/2, guessesY(jj), guessesX(jj), halfBoxWidth, sz);
    end
    
    % Loop for all good fits in the current frame
    for jj = 1:numel(r2)
        img = makeBox(img, 1, gfY(jj),gfX(jj), halfBoxWidth,sz);
    end
    
    % replicate grayscale to create RGB color images
    img = img(:,:,ones(1,3));
    
    % Loop for all tracks in the current frame
    for jj = 1:numel(r3)
        img = makeBox(img,cm(mod(trackFile(r3(jj),1),size(cm,1)-1)+1,:),...
            trY(jj),trX(jj),halfBoxWidth,sz);
    end
    
    %     imshow(img)
    %     v=getframe(gca);
    %     imwrite(currentframe2, [vidloc(1:end-4),'_Viewfits.tif'], 'writemode', 'append');
    writeVideo(wobj,img)
    %     writeVideo(wobj,v);
end
% close(wobj);
end

function currentframe=makeBox(currentframe,curr_fr_maxint,...
    p1,p2,half_symbol_size,sz)
b = p2-half_symbol_size;
t = p2+half_symbol_size;
l = p1-half_symbol_size;
r = p1+half_symbol_size;

b(b<1)=1; l(l<1)=1;
t(t>sz(1))=sz(1); r(r>sz(2))=sz(2);

if isnan(p1)||isempty(p1)||b>sz(1)||l>sz(2)||t<1||r<1
    return
end

for ii=1:size(currentframe,3)
    currentframe(b:t,[l,r],ii)=curr_fr_maxint(ii);
    currentframe([b,t],l:r,ii)=curr_fr_maxint(ii);
end
end

function phaseMask=selectCells(phaseMask,img)
% Let user click on the phase mask of cell images to decide which cell to
% analyze subsequently

%#ok<*AGROW>

h=figure;
c=onCleanup(@()close(h));
vSize = size(phaseMask);

subplot(121);
imshow(img,[]);
subplot(122);
imshow(phaseMask~=0,[]);

hold all

% highlight the outlines of the cells
rProp=regionprops(phaseMask,'Convexhull');
for i=1:numel(rProp)
    plot(rProp(i,1).ConvexHull(:,1),rProp(i,1).ConvexHull(:,2),'c-','linewidth',2)
end

set(gcf,'NextPlot','add');
axes
h1 = title('Click the cells to be analyzed. Enter to Proceed.');
set(gca,'Visible','off');
set(h1,'Visible','on');

% Let the user pick cells with good shapes
keyInput=0; clicksX=[]; clicksY=[];
while keyInput~=121
    [clickY,clickX,keyInput]=ginput(1);
    plot(clickY,clickX,'m*','markersize',12);
    clicksX=cat(1,clicksX,clickX);
    clicksY=cat(1,clicksY,clickY);
end

% If the user did click on something
if ~isempty(clicksX)
    pmID = findPmID(phaseMask,vSize,cat(2,clicksX,clicksY));
    phaseMask(~ismember(phaseMask,pmID)) = 0;
end
end

function pmID = findPmID(phaseMask,imSize,points)
pmID = zeros(size(points,1),1);

if numel(pmID)>0
    p1 = points(:,1) + 1 < imSize(1);
    p2 = points(:,2) + 1 < imSize(2);
    p3 = points(:,1)>0|points(:,2)>0;
    
    g1 = p1&p2&p3;
    
    points = ceil(points);
    
    pmID(g1) = phaseMask(sub2ind(imSize,points(g1,1),points(g1,2)));
end
end