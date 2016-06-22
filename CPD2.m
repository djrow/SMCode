function out = CPD2(mainfold,varargin)
% tracking data in the analysis files must be in pixels. the legacy
% "trackfile" variable has columns 4 and 5 in pixels.
%
% this histograms the bootstrapped results for the second diffusion
% coefficient (if it exists for the chosen fitting model):
% hist(out.fittedParameters_B(:,out.dID(1)))
%
% dID and aID give the indices for the diffusion coefficeints
% and population amplitudes in the fitted parameters vector.

opts = optimset('Display','off');
%#ok<*PFBNS>

nY = @(m,v,x)exp(-(x-m).^2/2/v)/sqrt(2*pi*v);

if ischar(mainfold)
    %% find the relevant files
    display('Select the analysis files.')
    [datalist,dataloc,~]=uigetfile([mainfold, filesep, '*.mat'],...
        'Matlab Files','multiselect','on');
    
    if ~dataloc
        display('no data selected')
        return
    end
    if ~iscell(datalist)
        datalist = {datalist};
    end
    
    tr=cell(1,numel(datalist)); tic
    for ii=1:numel(datalist)
        m=matfile(fullfile(dataloc,datalist{ii}));
        tr{ii}=m.trackfile;
    end; toc
else
    tr=mainfold;
    if ~iscell(tr)
        tr={tr};
    end
end

%% Default analysis parameters
anProp.nMobile = 4;     % number of diffusive populations
anProp.immBool = 0;     % presence or absence of stationary population
anProp.tFrame = .04;    % camera integration time in seconds
anProp.pixSize = .049;  % pixel size in microns
anProp.minTau = 1;      % minimum time lag in frames
anProp.maxTau = 5;      % maximum time lag in frames
anProp.stepSizeLP = 1;  % percentage of step sizes used for each tau
anProp.confBool = 0;    % confined or unconfined diffusion
anProp.globBool = 1;    % global or local cpd fit
anProp.overBool = 0;    % use overlapping or non-overlapping displacements?
anProp.plotBool = 1;    % plot output or not
anProp.dim = 2;         % 1d or 2d diffusion analysis
anProp.whichDim = 1;    % for 1d diffusion, which dimension to consider
anProp.rotAngle = pi/3; % for 1d diffusion, clockwise angle to rotate the trajectory
anProp.bootNum = 100;   % number of bootstraps

fNames=fieldnames(anProp);

% if any analysis parameters are included as inputs, change the analysis
% parameters mentioned
if nargin>1&&rem(nargin,2)==1
    for ii=1:2:nargin-1
        whichField = strcmp(fNames,varargin{ii});
        
        if all(~whichField)
            warning(['Check spelling. ', ...
                'Parameter change may have not occurred'])
        end
        
        eval([fNames{whichField} ' = varargin{ii+1};'])
        eval(['anProp.' fNames{whichField} ' = ' fNames{whichField},';'])
    end
elseif nargin>1
    warning('use paired inputs')
    dOut = [];
    return
end

% partial 2d cpd function
c2=@(x,y,p)p*exp(-x./y);

% partial 1d cpd function
c1=@(x,y,p)p*erf(sqrt(x./2/y));

% 2d confined msd function
m2=@(t,p)4*p(1)*t+p(2);

% 1d unconfined msd function
m1=@(t,p)2*p(1)*t+p(2);

%% calculate and accumulate squared step sizes
for kk=1:numel(tr)
    trackNums = unique(tr{kk}(:,1))';
    
    if ~isempty(trackNums)
        for ii = trackNums
            tracks = tr{kk}(tr{kk}(:,1)==ii,[2,4,5]);
            
            % fill in the time holes with nans
            fixedTrack = nan(max(tracks(:,1)),size(tracks,2));
            fixedTrack(tracks(:,1),:) = tracks;
            
            % remove leading nans
            fixedTrack(1:find(all(isnan(fixedTrack),2)==0,1,'first')-1,:) = [];
            
            nLocs = size(fixedTrack,1);
            for jj=1:anProp.maxTau
                if anProp.overBool
                    indvec1=jj+1:nLocs;
                    indvec2=1:nLocs-jj;
                elseif ~anProp.overBool
                    indvec2=1:jj:nLocs;
                    indvec1=indvec2(1:end-1);
                    indvec2=indvec2(2:end);
                end
                
                % calculate squared step sizes
                if anProp.dim == 2
                    allSqSteps{kk,ii,jj}=nansum( (fixedTrack(indvec1,[2,3]) - ...
                        fixedTrack(indvec2,[2,3])).^2, 2);
                elseif anProp.dim == 1
                    allSqSteps{kk,ii,jj}=(fixedTrack(indvec1,[2,3]) - ...
                        fixedTrack(indvec2,[2,3])).^2;
                end
            end
        end
    end
end

if anProp.dim == 2
    sqSteps=cell(anProp.maxTau,1);
    for ii=1:anProp.maxTau
        wSteps = cat(1,allSqSteps{:,:,ii});
        sqSteps{ii}=sort(wSteps(wSteps > eps));     % nansum puts zeros where there were nans
    end
    oRanks=cellfun(@(x)linspace(0,1,numel(x))',sqSteps,'uniformoutput',0);
else
    rMat = [cos(anProp.th),-sin(anProp.th);sin(anProp.th),cos(anProp.th)];
    
    sqSteps=cell(anProp.maxTau,1);
    for ii=1:anProp.maxTau
        y=sort(cat(1,allSqSteps{:,ii}));
        rPoints=y*rMat;
        sqSteps{ii}=rPoints(:,anProp.whichDim);
    end
    
    oRanks=cellfun(@(x)linspace(0,1,size(x,1))',sqSteps,'uniformoutput',0);
end
nSteps = cellfun(@numel,sqSteps,'uniformoutput',0);
nTotal = sum(cat(1,nSteps{:}));

%% fitting function selection
funFinds = cpdFunFinder2(anProp)
cpdFun = funFinds.cpdFun;
msdFun = funFinds.msdFun;
pStart = funFinds.pStart;
bounds = funFinds.bounds;
dID = funFinds.dID;
aID = funFinds.aID;

% time lag domain in seconds
tau = (1:anProp.maxTau)'*anProp.tFrame;
minTau = anProp.minTau;
if anProp.globBool
    %% GLOBAL FITTING
    linCell=@(x)cat(1,x{:});
    fHandle=@(p,tau,sqSteps,ranks)linCell(...
        cellfun(@(x,y)x-y,...
        cellfun(@(x,y)cpdFun(x,y,p),...
        sqSteps,num2cell(msdFun(tau,p),2),'uniformoutput',0),...
        ranks,'uniformoutput',0));
    eHandle=@(p,tau,sqSteps)cellfun(@(x,y)cpdFun(x,y,p),...
        sqSteps,num2cell(msdFun(tau,p),2),'uniformoutput',0);
    
    % fit to original data
    y=cellfun(@(x)x*anProp.pixSize.^2,sqSteps,'uniformoutput',0);
    y=cellfun(@(x,y)x(1:floor(y*anProp.stepSizeLP)),y,nSteps,'uniformoutput',0);
    r=cellfun(@(x,y)x(1:floor(y*anProp.stepSizeLP)),oRanks,nSteps,'uniformoutput',0);
    
    [fP_nB,~,r_nB] = lsqnonlin(@(p)fHandle(...
        p,tau(minTau:end),y(minTau:end),r(minTau:end)),...
        pStart{1},bounds{1},bounds{2},opts); disp(fP_nB)
    
    % bootstrap
    fP_B = zeros(anProp.bootNum,numel(pStart{1}));
    %     r_B = zeros(nTotal,anProp.bootNum);
    parfor kk = 1:anProp.bootNum
        % converts track positions from pixels to microns and resamples at the same
        % time
        y1=cellfun(@(x,y)sort(x(randsample(y,y,1))*anProp.pixSize.^2),sqSteps,nSteps,'uniformoutput',0);
        y1=cellfun(@(x,y)x(1:floor(y*anProp.stepSizeLP)),y1,nSteps,'uniformoutput',0);
        r1=cellfun(@(x,y)x(1:floor(y*anProp.stepSizeLP)),oRanks,nSteps,'uniformoutput',0);
        
        fP_B(kk,:) = lsqnonlin(@(p)fHandle(...
            p,tau(minTau:end),y1(minTau:end),r1(minTau:end)),...
            pStart{1},bounds{1},bounds{2},opts);
    end
    disp(fP_B(1,:))
    BIC = nTotal*log(mean(r_nB.^2,1)) + numel(pStart{1})*log(nTotal);
    
elseif ~anProp.globBool
    %% LOCAL FITTING
    cpdLB = bounds{1};
    cpdUB = bounds{2};
    msdLB = bounds{3};
    msdUB = bounds{4};
    cpdStart = pStart{1}; cpdStart = cpdStart(ones(1,anProp.maxTau),:);
    msdStart = pStart{2};
    
    % fit the cpds (original data)
    y = cellfun(@(x)x*anProp.pixSize.^2,sqSteps,'uniformoutput',0);
    r = [];
    for mm = 1:anProp.maxTau
        [fP1_nB(:,mm),~,rTemp] = lsqcurvefit(@(p,x)cpdFun(x,p),...
            cpdStart(mm,:),y{mm},oRanks{mm},cpdLB,cpdUB,opts);
        r = cat(1,r,rTemp);
    end
    r_nB = r;
    
    % fit the msd curves
    for mm = 1 : anProp.nMobile
        fP2_nB(:,mm) = lsqnonlin(@(p)(msdFun(tau,p)-fP1_nB((mm-1)*2+1,:)') .* [nSteps{:}]'/nTotal,...
            msdStart,msdLB,msdUB,opts);
    end
    
    %     for ii = 1:anProp.maxTau
    %         for jj = 1:anProp.nMobile
    %             cpdStart(ii,jj) = mean(y{ii})/10^(jj-1);
    %         end
    %     end
    
    fP2_B = zeros(anProp.bootNum,numel(msdStart),anProp.nMobile);
    parfor kk = 1:anProp.bootNum
        y=cellfun(@(x,y)sort(x(randsample(y,y,1))*anProp.pixSize.^2),sqSteps,nSteps,...
            'uniformoutput',0);
        
        % fit the cpd curves to get the msd values
        fParams1 = zeros(numel(cpdStart(1,:)),anProp.maxTau);
        for mm = 1:anProp.maxTau
            fParams1(:,mm) = lsqcurvefit(@(p,x)cpdFun(x,p),...
                cpdStart(mm,:),y{mm},oRanks{mm},cpdLB,cpdUB,opts);
        end
        
        fP1_B(:,:,kk) = fParams1;
        
        % fit msds to get diffusion coefficients
        temp = zeros(numel(msdStart),anProp.nMobile);
        for mm = 1 : anProp.nMobile
            % unweighted fit
            %         fParams2(:,mm)=lsqcurvefit(@(p,x)msdFun(x,p),msdStart,tau,fParams1(mm,:)',...
            %             msdLB,msdUB,opts);
            % weighted fit
            temp(:,mm) = lsqnonlin(@(p)(msdFun(tau,p)-fParams1(mm,:)') .* [nSteps{:}]'/nTotal,...
                msdStart,msdLB,msdUB,opts);
        end
        fP2_B(kk,:,:) = temp;
    end
    
    % calculate BIC
    BIC = nTotal * log(mean(r_nB.^2,1)) + ...
        anProp.maxTau * numel(pStart{1}) * log(nTotal);
end

out.BIC = BIC;
if anProp.globBool
    out.fittedParameters_B = fP_B;
    out.fittedParameters_nB = fP_nB;
    out.combinedResiduals_nB = r_nB;
    out.dID = dID;
    out.aID = aID;
else
    out.cpdParameters_B = fP1_B;
    out.cpdParameters_nB = fP1_nB;
    out.msdParameters_B = fP2_B;
    out.msdParameters_nB = fP2_nB;
    out.cpdResiduals_nB = r_nB;
end
end