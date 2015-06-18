function [dmeas,dmeas95]=CPD(trfile)
% computes msds from cumulative probability step size distributions, and
% then fits the msds to a msd model to estimate the diffusion coefficients
% 
% to use, run this without inputs and select the analysis files or run this
% with one input, which is the particular tracking file you want to
% analyze.

max_tau=5;          % maximum time lag in frames
inttime=.04;        % integration time in seconds

% number and type of terms in the CPD fit
nMobile=2;
nImmobile=0;
cpdstart=[];        % leave blank unless you know what you're doing

% mean squared displacement model for MSD fit
msdModel='2d square confinement';
msdstart=[];        % leave blank unless you know what you're doing

% use overlapping displacements?
yesoverlap=0;

% minimum track length
minTrLength=5;

% 2 or 1 dimensional diffusion
dim=2;

% choose the functions for fitting
[cpdfun,cpdstart]=cpd_function(nMobile,nImmobile,cpdstart);
[msdfun,msdstart]=msd_function(msdModel,msdstart);

% which taus to plot in the cpd fit result figure
plot_tau=[1,2,3];

% get the locations and names of all the analysis files
if ~nargin
    [trFileName,trFileLoc,idf]=uigetfile({'*.mat'},'Select analysis files.',...
        'MultiSelect','on');
    if ~idf
        display('no files chosen.')
        return
    end
    % convert the name to a cell array if only one file is chosen
    if ~iscell(trFileName)
        trFileName={trFileName};
    end
else
    % if there's an input, don't load anything
    trFileName={'manual input'};
end

temp=cell(1,max_tau); counter=0;
for kk=1:numel(trFileName)
    
    % if there's no input, load the files, if there's an input, just use
    % the input as the tracking file
    if ~nargin
        m1=matfile([trFileLoc,filesep,trFileName{kk}]);
        try
            trfile=m1.trfile;
        catch
            display([trFileName{kk} ' does not include tracking data. Skipping'])
            continue
        end
    end
    
    % loop through tracks
    for ii=unique(trfile(:,1))'
        counter=counter+1;
        
        % select current track from the array of all tracks
        trackii=trfile(trfile(:,1)==ii,:);
        
        % fill in the time holes with nans. these two lines are genius!
        fixedtrack=nan(max(trackii(:,2)),size(trackii,2));
        fixedtrack(trackii(:,2),:)=trackii;
        
        for jj=1:max_tau
            if yesoverlap                       % overlapping frame pairs
                indvec1=jj:size(fixedtrack,1);
                indvec2=1:size(fixedtrack,1)-jj;
            else                                % nonoverlapping frame pairs
                indvec2=1:jj:size(fixedtrack,1);
                indvec1=indvec2(1:end-1);
                indvec2=indvec2(2:end);
            end
            
            % nansum because there are nans as placeholders
            temp{counter,jj}=nansum((fixedtrack(indvec1,4:5)-...
                fixedtrack(indvec2,4:5)).^2,2);
        end
    end
end

sqsteps=cell(max_tau,1);
for ii=1:max_tau
    sqsteps{ii}=sort(cat(1,temp{:,ii}));
end
ranks=cellfun(@(x)linspace(0,1,numel(x))',sqsteps,'uniformoutput',0);

color_ind=0;
for ii=1:max_tau
    if numel(ranks{ii})<minTrLength
        continue
    end
    
    mdl=fitnlm(sqsteps{ii},ranks{ii},cpdfun,cpdstart);
    msds(ii,:)=mdl.Coefficients{:,1};
    
    cpdfxn=cpdfun(msds(ii),sqsteps{ii});
    residual=ranks{ii}-cpdfxn(:);
    
    % Create colormap for plotting raw CPD data
    cmap=hsv(numel(plot_tau));
    
    %  plotting CPD at specified tau
    if ismember(ii,plot_tau)
        color_ind=color_ind+1;
        
        % Plot CPDs
        subplot(50,1,1:40)
        semilogx(sqsteps{ii},ranks{ii},'.','Color',cmap(color_ind,:),...
            'MarkerSize',8);
        hold all
        
        % Plot fitted lines
        semilogx(sqsteps{ii},cpdfxn,'-','color','k','Linewidth',2,...
            'HandleVisibility','off')
        
        % Plot residuals
        subplot(50,1,41:50)
        semilogx(sqsteps{ii},residual,'.','Color',cmap(color_ind,:),...
            'MarkerSize',4);
        hold all
        
        subplot(50,1,1:40)
        set(gca,'XTickLabel',[])
    end
end

tau=inttime:inttime:maxtau*inttime;
dmeas=nan(nMobile,numel(msdstart));
for ii=nMobile+nImmobile:numel(cpdstart)-nImmobile
    y=msds(:,ii);
    if ~sum(~isnan(y))
        display(['diffusive population number ', ...
            num2str(ii-nMobile-nImmobil+1), ' has no msd values'])
        dmeas(ii,:)=nan(size(msdstart));
        
        continue
    end
    mdl=fitnlm(tau(~isnan(y)),y(~isnan(y)),msdfun,msdstart);
    dmeas(ii,:)=mdl.Coefficients{:,2};
    dmeas95(ii,:)=diff(coefCI(mdl),[],2);
end
dmeas=abs(dmeas);

end

%% ------------------------------------------------------------------------
%  Available CPD Model Functions
%  ------------------------------------------------------------------------
function [fhandle,pstart]=cpd_function(nMobile,nImmobile,pstart)
switch nMobile
    case 1
        switch nImmobile
            case 0
                fhandle=@(p,sqr)1-...
                    exp(-sqr/abs(p(1)));
                if isempty(pstart)
                    pstart=1;
                end
            case 1
                fhandle=@(p,sqr)1-...
                    abs(p(1))*exp(-sqr/(abs(p(2))+abs(p(3))))-...
                    (1-abs(p(1)))*exp(-sqr/abs(p(3)));
                if isempty(pstart)
                    pstart=[.5,1,eps];
                end
        end
    case 2
        switch nImmobile
            case 0
                fhandle=@(p,sqr)1-...
                    abs(p(1))*exp(-sqr/abs(p(2)))-...
                    (1-abs(p(1)))*exp(-sqr/abs(p(3)));
                if isempty(pstart)
                    pstart=[.5,1,.1];
                end
            case 1
                fhandle=@(p,sqr)1-...
                    abs(p(1))*exp(-sqr/(abs(p(3))+abs(p(5))))-...
                    abs(p(2))*exp(-sqr/(abs(p(4))+abs(p(5))))-...
                    (1-abs(p(1))-abs(p(2)))*exp(-sqr/(abs(p(5))));
                if isempty(pstart)
                    pstart=[.3,.3,1,.1,eps];
                end
        end
    case 3
        switch nImmobile
            case 0
                fhandle=@(p,sqr)1-...
                    abs(p(1))*exp(-sqr/abs(p(3)))-...
                    abs(p(2))*exp(-sqr/abs(p(4)))-...
                    (1-abs(p(1))-abs(p(2)))*exp(-sqr/abs(p(5)));
                if isempty(pstart)
                    pstart=[.3,.3,1,.1,.01];
                end
            case 1
                fhandle=@(p,sqr)1-...
                    abs(p(1))*exp(-sqr/(abs(p(4))+abs(p(7))))-...
                    abs(p(2))*exp(-sqr/(abs(p(5))+abs(p(7))))-...
                    abs(p(3))*exp(-sqr/(abs(p(6))+abs(p(7))))-...
                    (1-abs(p(1))-abs(p(2))-abs(p(3)))*exp(-sqr/(abs(p(7))));
                if isempty(pstart)
                    pstart=[.25,.25,.25,1,.1,.01,eps];
                end
        end
end
end

%% ------------------------------------------------------------------------
%  Available MSD Model Functions
%  ------------------------------------------------------------------------
function [fhandle,pstart]=msd_function(model,pstart)

switch model
    case '1d linear'
        fhandle=@(p,x) 2*abs(p(1))*x;
        if isempty(pstart)
            pstart=.1;
        end
    case '2d linear'
        fhandle=@(p,x) 4*abs(p(1))*x;
        if isempty(pstart)
            pstart=.1;
        end
    case '1d square confinement'
        fhandle=@(p,x) longmsd(p,x)/2;
        if isempty(pstart)
            pstart=[1,.1,0];
        end
    case '2d square confinement'
        fhandle=@longmsd;
        if isempty(pstart)
            pstart=[1,.1,0];
        end
    otherwise
        display('mismelled the confinement model type')
        fhandle=[];
        pstart=[];
end
end

%% square confinement model
function z=longmsd(p,x)
% global camerasd
p=abs(p);

l=p(1); d=p(2); tau=x;

summedterm=@(t,d,l,n)1/n^4*exp(-(n*pi/l).^2*d*t);

temp=eps*ones(size(tau));% counter=0;
for ii=1:2:2*400-1
    s=summedterm(tau,d,l,ii);
    if sum(s./temp)<1e-10
        break
    end
    temp=temp+s;
    %     counter=counter+1;
end
% display(counter)
z=l^2/3*(1-96/pi^4*temp)+p(3); %.0384;
end