function [dOut] = CPD2(mainfold,varargin)
% tracking data in the analysis files must be in pixels.

opts = optimset('Display','off');

if ischar(mainfold)
    %% find the relevant files
    display('Select the analysis files.')
    [datalist,dataloc,~]=uigetfile([mainfold, filesep, '*.mat'],...
        'Matlab Files','multiselect','on');
    
    if ~dataloc
        display('no data selected')
        dOut = [];
        return
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
anProp.nMobile = 1;     % number of diffusive populations
anProp.immBool = 0;     % presence or absence of stationary population
anProp.tFrame = .04;    % camera integration time in seconds
anProp.pixSize = .049;  % pixel size in microns
anProp.maxTau = 5;      % maximum time lag in frames
anProp.confBool = 0;    % confined or unconfined diffusion
anProp.globBool = 1;    % global or local cpd fit
anProp.overBool = 1;    % use overlapping or non-overlapping displacements?
anProp.plotBool = 0;    % plot output or not
anProp.dim = 2;         % 1d or 2d diffusion analysis
anProp.whichDim = 1;    % for 1d diffusion, which dimension to consider
anProp.rotAngle = pi/3; % for 1d diffusion, clockwise angle to rotate the trajectory
anProp.bootNum = 200;   % number of bootstrap iterations

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
                    allSqSteps{ii,jj}=nansum( (fixedTrack(indvec1,[2,3]) - ...
                        fixedTrack(indvec2,[2,3])).^2, 2);
                elseif anProp.dim == 1
                    allSqSteps{ii,jj}=(fixedTrack(indvec1,[2,3]) - ...
                        fixedTrack(indvec2,[2,3])).^2;
                end
            end
        end
    end
end

if anProp.dim == 2
    sqSteps=cell(anProp.maxTau,1);
    for ii=1:anProp.maxTau
        sqSteps{ii}=sort(cat(1,allSqSteps{:,ii}));
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

%% fitting function selection
[cpdFun,msdFun,pStart,bounds,dID] = cpdFunFinder(anProp);

% time lag domain in seconds
tau = (1:anProp.maxTau)'*anProp.tFrame;

dOut = zeros(numel(dID),anProp.bootNum);
h=waitbar(0,'boot reps');
c=onCleanup(@()close(h));
for kk = 1:anProp.bootNum
    waitbar(kk/anProp.bootNum)
    
    % converts track positions from pixels to microns
    y=cellfun(@(x,y)sort(x(randsample(y,y,1))*anProp.pixSize.^2),sqSteps,nSteps,'uniformoutput',0);
    %         y=cellfun(@(x)x*anProp.pixSize.^2,sqSteps,'uniformoutput',0);
    
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
        
        fParams=lsqnonlin(@(p)fHandle(p,tau,y,oRanks),...
            pStart{1},bounds{1},bounds{2},opts);
        
        temp=cellfun(@(x)x*anProp.pixSize.^2,sqSteps,'uniformoutput',0);
        resids = cellfun(@(x,y)x-y,oRanks,eHandle(fParams,tau,temp),'uniformoutput',0);
        
        dOut(:,kk)=fParams(dID);
    elseif ~anProp.globBool
        %% LOCAL FITTING
        % fit the cpd curves to get the msd values
        for ii=1:anProp.maxTau
            cpdLB = bounds{1};
            cpdUB = bounds{2};
            cpdStart = pStart{1};
            
            % this loop replaces nans
            for jj = 1:anProp.nMobile
                cpdStart(jj) = mean(y{ii})/(10^(jj-1))*anProp.tFrame;
            end
            
            msds(:,ii)=lsqcurvefit(@(p,x)cpdFun(x,p),...
                cpdStart,y{ii},oRanks{ii},cpdLB,cpdUB,opts);
            resids{1,ii} = cpdFun(y{ii},msds(:,ii)) - oRanks{ii};
        end
        
        % fit msds to get diffusion coefficient
        for ii = 1 : anProp.nMobile
            y=msds(ii,:)';
            msdStart = pStart{2};
            msdLB = bounds{3};
            msdUB = bounds{4};
            
            fParams(ii,:)=lsqcurvefit(@(p,x)msdFun(x,p),msdStart,tau,y,...
                msdLB,msdUB,opts);
            resids{2,ii} = msdFun(tau,fParams(ii,:)) - y;
            
            dOut(ii,kk)=fParams(ii,1);
        end        
    end
end

%% plot results
if ~anProp.globBool&&anProp.plotBool
    a=figure;
    b=figure;
    
    figure(a)
    for ii = 1:anProp.maxTau
        steps = sqSteps{ii}*anProp.pixSize.^2;
        ranks = oRanks{ii};
        
        % data
        ind=sub2ind([anProp.maxTau,2],ii,1);
        subplot(2,anProp.maxTau,ind);
        title(['data for tau = ' num2str(ii)])
        
        % fittign result
        loglog(steps,cat(2,ranks,cpdFun(steps,msds(:,ii))),'.')
        
        % residuals
        ind=sub2ind([anProp.maxTau,2],ii,2);
        subplot(2,anProp.maxTau,ind); title('residuals')
        loglog(steps,abs(resids{1,ii}),'.')
    end
    
    figure(b)
    for ii = 1:anProp.nMobile
        % data
        subplot(2,anProp.nMobile,(ii-1)*2+1)
        plot(tau,msds(ii,:),'o'); title('data'); hold on
        plot(tau,msdFun(tau,fParams(ii,:))); hold off
        
        % residuals
        subplot(2,anProp.nMobile,(ii-1)*2+2);
        plot(tau,resids{2,ii}); title('residuals')
    end
elseif anProp.plotBool
    figure
    y = cellfun(@times,sqSteps,num2cell(anProp.pixSize.^2*ones(anProp.maxTau,1)),'uniformoutput',0);
    fRes = eHandle(fParams,tau,y);
    
    axLims{1}(1:2) = [min(cellfun(@(x)min(x),y)),...
        max(cellfun(@(x)max(x),y))];
    axLims{1}(3:4) = [min(cellfun(@(x)min(x),oRanks)),max(cellfun(@(x)max(x),oRanks))];
    
    axLims{2}(1:2) = [min(cellfun(@(x)min(abs(x)),resids)),...
        max(cellfun(@(x)max(abs(x)),resids))];
    for ii = 1:anProp.maxTau
        steps = sqSteps{ii}*anProp.pixSize.^2;
        ranks = oRanks{ii};
        
        % data
        ind=sub2ind([anProp.maxTau,2],ii,1);
        subplot(2,anProp.maxTau,ind);
        loglog(steps,cat(2,ranks,fRes{ii}),'.')
        title(['data for tau = ' num2str(ii)])
        axis(axLims{1})
        
        % residuals
        ind=sub2ind([anProp.maxTau,2],ii,2);
        subplot(2,anProp.maxTau,ind)
        loglog(steps,abs(resids{ii}),'.');
        title('residuals')
        axis([axLims{1}(1:2),axLims{2}])
    end
    hold off
end

figure
for ii=1:numel(dID)
%     subplot(1,numel(dID),ii)
    ecdf(dOut(ii,:)); hold all
end; hold off
xlabel('d in microns^2/s')
end