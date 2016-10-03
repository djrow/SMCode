function [out,sqSteps] = CPD2(mainfold,varargin)
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
    
    tr=cell(1,numel(datalist));
    for ii=1:numel(datalist)
        m=matfile(fullfile(dataloc,datalist{ii}));
        tr{ii}=m.trackfile;
    end;
else
    tr=mainfold;
    if ~iscell(tr)
        tr={tr};
    end
end

%% Default analysis parameters
anProp.nMobile = 1;     % number of diffusive populations
anProp.immBool = 0;     % presence or absence of stationary population
anProp.tFrame = .0104;    % camera integration time in seconds
anProp.pixSize = .049;  % pixel size in microns
anProp.minTau = 1;      % minimum time lag in frames
anProp.maxTau = 5;      % maximum time lag in frames
anProp.stepSizeLP = 1;  % percentage of step sizes used for each tau
anProp.confBool = 0;    % confined or unconfined diffusion
anProp.globBool = 1;    % global or local cpd fit
anProp.overBool = 0;    % use overlapping or non-overlapping displacements?
anProp.plotBool = 0;    % plot output or not
anProp.dim = 2;         % 1d or 2d diffusion analysis
anProp.whichDim = 1;    % for 1d diffusion, which dimension to consider
anProp.rotAngle = 0;    % for 1d diffusion, clockwise angle to rotate the trajectory
anProp.bootNum = 300;   % number of bootstraps

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
    Out = [];
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
    if isempty(tr{kk})
        continue
    end
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
    rMat = [cos(anProp.rotAngle),-sin(anProp.rotAngle);...
        sin(anProp.rotAngle),cos(anProp.rotAngle)];
    
    sqSteps=cell(anProp.maxTau,1);
    for ii=1:anProp.maxTau
        y=sort(cat(1,allSqSteps{:,ii}));
        rPoints=y*rMat;
        sqSteps{ii}=rPoints(:,anProp.whichDim);
    end
    
    oRanks=cellfun(@(x)linspace(0,1,size(x,1))',sqSteps,'uniformoutput',0);
end
sqSteps = cellfun(@(x)x*anProp.pixSize.^2,sqSteps,'uniformoutput',0);
nSteps = cellfun(@numel,sqSteps,'uniformoutput',0);
nTotal = sum(cat(1,nSteps{:}));

%% fitting function selection
funFinds = cpdFunFinder2(anProp);
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
    y=sqSteps;
    y=cellfun(@(x,y)x(1:floor(y*anProp.stepSizeLP)),y,nSteps,'uniformoutput',0);
    r=cellfun(@(x,y)x(1:floor(y*anProp.stepSizeLP)),oRanks,nSteps,'uniformoutput',0);
    
    % evenly space the diffusion coefficient starting values
    if numel(dID)>1
        pStart{1}(dID) = logspace(-5,0,numel(dID));
        %         pStart{1}(dID) = [1e-4,1e-1];
    else
        pStart{1}(dID) = .01;
    end
    
    % uninformed amplitude guesses
    pStart{1}(aID) = 1/(numel(aID)+1);
    
    [fP_nB,~,r_nB] = lsqnonlin(@(p)fHandle(...
        p,tau(minTau:end),y(minTau:end),r(minTau:end)),...
        pStart{1},bounds{1},bounds{2},opts);
    
    rTemp = r_nB;
    for ii = anProp.minTau:anProp.maxTau
        residCell{ii} = rTemp(1:nSteps{ii});
        rTemp(1:nSteps{ii}) = [];
    end
    
    % bootstrap
    fP_B = zeros(anProp.bootNum,numel(pStart{1}));
    ci = zeros(anProp.bootNum,numel(pStart{1}));
    %     r_B = zeros(nTotal,anProp.bootNum);
    parfor kk = 1:anProp.bootNum
        % resamples with replacement and selects only the steps smaller than the
        % threshold
        y1=cellfun(@(x,y)sort(x(randsample(y,y,1))),sqSteps,nSteps,'uniformoutput',0);
        y1=cellfun(@(x,y)x(1:floor(y*anProp.stepSizeLP)),y1,nSteps,'uniformoutput',0);
        r1=cellfun(@(x,y)x(1:floor(y*anProp.stepSizeLP)),oRanks,nSteps,'uniformoutput',0);
        
        [fP_B(kk,:),~, ~,~,~,~, ~] = lsqnonlin(@(p)fHandle(...
            p,tau(minTau:end),y1(minTau:end),r1(minTau:end)),...
            pStart{1},bounds{1},bounds{2},opts);
        
        ci(kk,:) = zeros(numel(pStart{1}),1);
        %         ci(kk,:) = diff(nlparci(fP_B(kk,:),resid,'jacobian',jac),[],2);
    end
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
    y = sqSteps;

    if anProp.nMobile == 2
        cpdStart(:,1) = .02+(0:anProp.maxTau-1)*.036;
        cpdStart(:,2) = .008+(0:anProp.maxTau-1)*.014;
    end
    
    r_nB = [];
    for mm = anProp.minTau:anProp.maxTau
        [fP1_nB(:,mm),~,rTemp] = lsqcurvefit(@(p,x)cpdFun(x,p),...
            cpdStart(mm,:),y{mm},oRanks{mm},cpdLB,cpdUB,opts);
        r_nB = cat(1,r_nB,rTemp);
        residCell{mm} = rTemp;
    end
    
    % fit the msd curves
    for mm = 1 : anProp.nMobile
        %         fP2_nB(:,mm) = lsqcurvefit(@(p,x)msdFun(x,p),msdStart,msdLB,msdUB,opts);
        [fP2_nB(:,mm),~,~] = lsqnonlin(@(p)(...
            msdFun(tau(anProp.minTau:anProp.maxTau),p)-...
            fP1_nB(mm,anProp.minTau:anProp.maxTau)') .* [nSteps{anProp.minTau:anProp.maxTau}]'/nTotal,...
            msdStart,msdLB,msdUB,opts);
    end
    
    fP2_B = zeros(numel(msdStart),anProp.nMobile,anProp.bootNum);
    parfor kk = 1:anProp.bootNum
        y=cellfun(@(x,y)sort(x(randsample(y,y,1))),sqSteps,nSteps,...
            'uniformoutput',0);
        
        % fit the cpd curves to get the msd values
        fParams1 = zeros(numel(cpdStart(1,:)),anProp.maxTau);
        for mm = anProp.minTau:anProp.maxTau
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
            temp(:,mm) = lsqnonlin(@(p)(msdFun(tau(anProp.minTau:anProp.maxTau),p)-...
                fParams1(mm,anProp.minTau:anProp.maxTau)') .* ...
                [nSteps{anProp.minTau:anProp.maxTau}]'/nTotal,...
                msdStart,msdLB,msdUB,opts);
        end
        fP2_B(:,:,kk) = temp;
    end
    
    % calculate BIC
    BIC = nTotal * log(mean(r_nB.^2,1)) + ...
        anProp.maxTau * numel(pStart{1}) * log(nTotal);
end

out.BIC = BIC;

if anProp.nMobile > 1
    nPlots = 3;
else
    nPlots = 2;
end
if anProp.globBool
    try
        out.fittedParameters_B = fP_B;
    end
    out.fittedParameters_nB = fP_nB;
    out.cpdResiduals = residCell;
    out.sqSteps = sqSteps;
    out.dID = dID;
    out.aID = aID;
    out.ci = ci;
    if anProp.plotBool
        c=colormap('lines');
        c=c(1:7,:);
        
        subplot(nPlots,1,1)
        for ii = 1:numel(dID)
            h=histogram(fP_B(:,dID(ii)),'normalization','probability','displaystyle','stairs');
            set(h,'edgecolor',c(ii,:))
            hold on
        end
        title('Bootstrapped diffusion coefficients')
        set(gca,'xscale','log'); hold off
        
        subplot(nPlots,1,2)
        for ii = anProp.minTau:anProp.maxTau
            plot(sqSteps{ii},residCell{ii}); hold all
        end
        title('Original data CPD residuals');
        hold off
        
        if anProp.nMobile > 1
            subplot(nPlots,1,3)
            lastAmp = ones(anProp.bootNum,1);
            for ii = 1:numel(aID)
                h=histogram(fP_B(:,aID(ii)),'normalization','probability','displaystyle','stairs');
                set(h,'edgecolor',c(ii,:))
                hold on
                
                lastAmp = lastAmp-fP_B(:,aID(ii));
            end
            h=histogram(lastAmp,'normalization','probability','displaystyle','stairs');
            title('Bootstrapped population amplitudes')
            hold off
        end
    end
else
    try
        out.cpdParameters_B = fP1_B;
        out.msdParameters_B = fP2_B;
    end
    out.cpdParameters_nB = fP1_nB;
    out.msdParameters_nB = fP2_nB;
    out.cpdResiduals = residCell;
    out.sqSteps = sqSteps;
    
    if anProp.plotBool
        c=[0    0.4470    0.7410;
            0.8500    0.3250    0.0980;
            0.9290    0.6940    0.1250;
            0.4940    0.1840    0.5560];
        
        % diffusion coefficients
        subplot(nPlots,1,1)
        for ii = 1:size(fP2_B,2)
            h=histogram(squeeze(fP2_B(1,ii,:)),'normalization','probability','displaystyle','stairs');
            set(h,'edgecolor',c(ii,:))
            hold on
        end
        title('Bootstrapped diffusion coefficients')
        hold off
        
        % cpd residuals
        subplot(nPlots,1,2)
        for ii = anProp.minTau:anProp.maxTau
            plot(sqSteps{ii},residCell{ii}); hold all
        end
        set(gca,'xscale','log')
        title('Original data CPD residuals')
        hold off
        
        if anProp.nMobile > 1
            % population amplitudes
            subplot(nPlots,1,3)
            lastAmp = ones(anProp.bootNum,1);
            counter = 0;
            for ii = anProp.nMobile+1:size(fP1_B,1)
                counter = counter + 1;
                h=histogram(squeeze(mean(fP1_B(ii,:,:))),'normalization','probability','displaystyle','stairs');
                set(h,'edgecolor',c(counter,:))
                hold on
                
                lastAmp = lastAmp-squeeze(mean(fP1_B(ii,:,:)));
            end
            h=histogram(lastAmp,'normalization','probability','displaystyle','stairs');
            title('Bootstrapped mean population amplitudes')
            hold off
        end
    end
end

end