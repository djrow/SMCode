function x=CPD2(tracks,varargin)
% tracks is a cell array of tracks. each element i is [2 x ni]

%% Default analysis parameters
tFrame = .05;           % camera integration time in seconds
maxTau = 5;             % maximum time lag in frames
confBool = 1;           % confined or unconfined diffusion
overBool = 1;           % use overlapping or non-overlapping displacements?
dim = 2;                % 1d or 2d diffusion analysis
globBool = 1;           % global or local cpd fit
nMobile = 1;            % number of diffusive populations
immBool = 0;            % presence or absence of stationary population
plotBool = 1;           % plot output or not

% initialize analysis properties structure
anProp.nMobile = nMobile;
anProp.immBool = immBool;
anProp.tFrame = tFrame;
anProp.maxTau = maxTau;
anProp.confBool = confBool;
anProp.globBool = globBool;
anProp.overBool = overBool;
anProp.dim = dim;

% if any analysis parameters are included as inputs, change the analysis
% parameters mentioned
if nargin>1&&rem(nargin,2)==1
    fNames=fieldnames(anProp);
    for ii=2:2:nargin
        whichField = strcmp(fNames,varargin{ii-1});
        
        if all(~whichField)
            warning(['Check spelling. ', ...
                'Parameter change may have not occurred'])
        end
        
        eval([fNames{whichField} ' = varargin{ii}'])
        eval(['anProp.' fNames{whichField} ' = ' fNames{whichField},';'])
    end
elseif nargin>1
    warning('use paired inputs')
    dmeas = [];
end

%% calculate and accumulate squared step sizes
if ~iscell(tracks)
    tracks = {tracks};
end

allSqSteps=cell(1,maxTau);
for kk=1:numel(tracks)
    for jj=1:maxTau
        if overBool
            indvec1=jj+1:size(tracks{kk},2);
            indvec2=1:size(tracks{kk},2)-jj;
        elseif ~overBool
            indvec2=1:jj:size(tracks{kk},2);
            indvec1=indvec2(1:end-1);
            indvec2=indvec2(2:end);
        end
        
        if dim == 2
            allSqSteps{jj}=sum( (tracks{kk}(:,indvec1) - ...
                tracks{kk}(:,indvec2)).^2, 1);
        elseif dim == 1
            allSqSteps{jj}=(tracks{kk}(:,indvec1) - ...
                tracks{kk}(:,indvec2)).^2;
        end
    end
end

sqSteps=cell(maxTau,1);
for ii=1:maxTau
    sqSteps{ii}=sort(cat(1,allSqSteps{:,ii}));
end
oRanks=cellfun(@(x)linspace(0,1,numel(x)),sqSteps,'uniformoutput',0);

%% fitting function selection
[cpdFun,msdFun,pStart,bounds]= ...
    cpdFunFinder(dim,nMobile,immBool,confBool,globBool);

% time lag domain
tau = (1:maxTau)'*tFrame;

if globBool
    linCell=@(x)cat(2,x{:});
    fHandle=@(p,tau,sqSteps,ranks)linCell(...
        cellfun(@(x,y)x-y,...
        cellfun(@(x,y)cpdFun(x,y,p),...
        sqSteps,num2cell(msdFun(tau,p),2),'uniformoutput',0),...
        ranks,'uniformoutput',0));
    eHandle=@(p,tau,sqSteps)cellfun(@(x,y)cpdFun(x,y,p),...
        sqSteps,num2cell(msdFun(tau,p),2),'uniformoutput',0);
    
    x=lsqnonlin(@(p)fHandle(p,tau,sqSteps,oRanks),...
        pStart{1},bounds{1},bounds{2});
    
    resids = cellfun(@(x,y)x-y,oRanks,eHandle(x,tau,sqSteps),'uniformoutput',0);
    
elseif ~globBool
    %% fit the cpd curves to get the msd values
    for ii=1:maxTau
        cpdLB = bounds{1};
        cpdUB = bounds{2};
        cpdStart = pStart{1};
        msds(:,ii)=lsqcurvefit(cpdFun,cpdStart,sqSteps{ii},oRanks{ii},...
            cpdLB,cpdUB);
        resids{1}(:,ii) = cpdFun(sqSteps{ii},msds(:,ii)) - oRanks{ii};
    end
    
    %% fit msds
    
    x=nan(nMobile,numel(pStart{2}));
    for ii = 1 : nMobile
        y=msds(ii,:);
        msdStart = pStart{2};
        msdLB = bounds{3};
        msdUB = bounds(4);
        
        x(ii,:)=lsqcurvefit(msdFun,msdStart,tau,y,...
            msdLB,msdUB,opts);
        resids{2}(:,ii) = msdFun(tau,x(ii,:)) - y;
    end
end

%% plot results
if ~globBool&&plotBool
    a=figure;
    b=figure;
    
    figure(a)
    for ii = 1:maxTau
        steps = sqSteps{ii};
        ranks = oRanks{ii};
        
        % data
        subplot(2,maxTau,(ii-1)*2+1); title('data');
        plot(steps,ranks)
        
        % fittine result
        plotyy(steps,cpdFun(steps,msds(:,ii)))
        
        % residuals
        subplot(2,maxTau,(ii-1)*2+2); title('residuals')
        plot(steps,resids{1}(:,ii))
    end
    
    figure(b)
    for ii = 1:nMobile
        % data
        subplot(2,nMobile,(ii-1)*2+1); title('data')
        plot(tau,msds(ii,:))
        
        % residuals
        subplot(2,nMobile,(ii-1)*2+2); title('residuals')
        plot(tau,resids{2}(:,ii))
    end
elseif plotBool
    fRes = eHandle(x,tau,sqSteps);
    
    axLims{1}(1:2) = [min(cellfun(@(x)min(x),oRanks)),max(cellfun(@(x)max(x),oRanks))];
    axLims{1}(3:4) = [min(cellfun(@(x)min(x),sqSteps)),max(cellfun(@(x)max(x),sqSteps))];
    
    axLims{2}(1:2) = [min(cellfun(@(x)min(abs(x)),resids)),...
        max(cellfun(@(x)max(abs(x)),resids))];
    for ii = 1:maxTau
        steps = sqSteps{ii};
        ranks = oRanks{ii};
        
        % data
        ind=sub2ind([maxTau,2],ii,1);
        subplot(2,maxTau,ind);
        loglog(steps,cat(1,ranks,fRes{ii}),'.')
        title(['data for tau = ' num2str(ii)])
        axis(axLims{1})
        
        % residuals
        ind=sub2ind([maxTau,2],ii,2);
        subplot(2,maxTau,ind); title('residuals')
        loglog(steps,abs(resids{ii}),'.')
        title('residuals')
        axis([axLims{1}(1:2),axLims{2}])
    end
    hold off
end
end