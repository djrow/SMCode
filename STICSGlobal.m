function outPut=STICSGlobal(v,nDiffs,immPop,diffModel)

% nDiffs=1;
% immPop=0;
% diffModel='confined';

% know what you're doing
pStart=[];

% filters to use when calculating the correlation functions
% 1:
% 2:
% 3: pad with mean
% 4: non-overlapping time lags;
% 5: pad with zeros
sticsFilters=[3,4];

% maximum frame separation to consider
maxTau=10;

% magnification of the microscope in micrometers/pixel.
mpp=.049;

% time from beginning of one frame to the beginning of the next in seconds
intTime=.04;

% rotation angle
th=pi/2;

% Choose the functions for fitting
if isempty(pStart)
    pStart=[.1,1,.01,.5,.01,...  % msd 1 parameters
        .01,.01,.01,...         % msd 2 parameters
        0,.01,0,0,.01,.01,.01];         % correlation function parameters
end
cpdLB([3,5,7,8,9,11,12])=-inf;
cpdLB([1,2,4,6,10,13,14,15])=0;
cpdUB=inf(1,50);

[corrFun,msdFun,pId]=sticsFunFinder2(nDiffs,immPop,diffModel,th);
pStart=pStart(pId);
cpdLB=cpdLB(pId);
cpdUB=cpdUB(pId);

% gooooooooooood luck figuring this out. nanLinCell linearizes a cell array
% into a single 1D vector. nanLinCell and longmsd both have to be aux
% functions in the code that fHandle and eHandle are used in.
fHandle=@(p,tau,corrDom,corrVal)nanLinCell(cellfun(@(x,y)x-y,...
    cellfun(@(x,y)corrFun(x,y,p),corrDom,num2cell(msdFun(tau,p),1),...
    'uniformoutput',0),corrVal,'uniformoutput',0));
eHandle=@(p,tau,corrDom)cellfun(@(x,y)corrFun(x,y,p),...
    corrDom,num2cell(msdFun(tau,p),1),'uniformoutput',0);

%% calculate correlations
display('correlating...')
tCorr=stics3(v,[],maxTau,sticsFilters);

%% fit correlations
wSizeP=size(tCorr);

% space correlations for each tau must be in a cell element by themselves
ind1=round(wSizeP(1)/2-wSizeP(1)/4):round(wSizeP(1)/2+wSizeP(1)/4);
ind2=round(wSizeP(2)/2-wSizeP(2)/4):round(wSizeP(2)/2+wSizeP(2)/4);
tCorr=squeeze(mat2cell(tCorr(ind1,ind2,2:end),...
    numel(ind1),numel(ind2),ones(1,maxTau)))';

% space domain
[X,Y]=ndgrid((ind1-mean(ind1))*mpp,(ind2-mean(ind2))*mpp);
corrDom=cat(3,X,Y);
corrDom={corrDom};
corrDom=corrDom(1,ones(1,maxTau));

% time domain
tau=(1:maxTau)*intTime;

% fit
% pStart(3)=max(tCorr{1}(:))-min(tCorr{1}(:));
display('fitting...')
fCoeff=lsqnonlin(@(p)fHandle(p,tau,corrDom,tCorr),pStart,cpdLB,cpdUB);

%% Plot the results
predCorr=eHandle(fCoeff,tau,corrDom);
startCorr=eHandle(pStart,tau,corrDom);
residCorr=cellfun(@(x,y)x-y,tCorr,predCorr,'uniformoutput',0);

pTau=3;
for kk=1:pTau
    subplot(3,pTau,kk); pcolor(tCorr{(kk-1)*3+1});
    axis image; shading flat; axis off
    cl=get(gca,'clim');
    
    subplot(3,pTau,kk+pTau); pcolor(residCorr{(kk-1)*3+1});
    axis image; shading flat; axis off
    cl=cl-min(cl)+nanmin(nanmin(residCorr{(kk-1)*3+1}));
    set(gca,'clim',cl);
    
    subplot(3,pTau,kk+2*pTau); pcolor(startCorr{(kk-1)*3+1});
    axis image; shading flat; axis off
    cl=cl-min(cl)+nanmin(nanmin(startCorr{(kk-1)*3+1}));
    set(gca,'clim',cl);
end

%% collect and label outputs
% non-existent fitting parameters exist as empty fields in outPut.
outPut.diffC1       =fCoeff(pId==1);    % diffusion coefficient 1
outPut.confX        =fCoeff(pId==2);    % confinement length 1
outPut.locNoise1    =fCoeff(pId==3);
outPut.confY        =fCoeff(pId==4);    % confinement length 2
outPut.locNoise2    =fCoeff(pId==5);
outPut.diffC2       =fCoeff(pId==6);    % diffusion coefficient 2
outPut.locNoise3    =fCoeff(pId==7);
outPut.locNoise4    =fCoeff(pId==8);
outPut.corrOffset   =fCoeff(pId==9);    % correlation offset
outPut.amp1         =fCoeff(pId==10);   % correlation amplitude 1
outPut.corrCenterX  =fCoeff(pId==11);   % correlation center 1
outPut.corrCenterY  =fCoeff(pId==12);   % correlation center 2
outPut.amp2         =fCoeff(pId==13);   % correlation amplitude 2
outPut.ampImm       =fCoeff(pId==14);   % correlation amplitude 3
outPut.immNoise     =fCoeff(pId==15);

end
function out=nanLinCell(in)
out=cat(1,in{:});
out=out(~isnan(out));
end