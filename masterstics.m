function output=masterstics(checkvals,mainfold)
% modes of use of this function

% 1. checkvals is a logical and mainfold is a string specifying the
% location of the folder containing the data

% 3.

output.msdp=[];

if ~exist('mainfold','var')
    mainfold=pwd;
end

% select cells from phase mask?
phasemasks=0;

% skip the selection of cells from the phase mask
skipselect=0;

% calculate correlation functions from data?
yesStics=0;

% fit the correlation function to estimate various observable phenomena
yesGauss=1;

% use local or global fit
localFit=0;

% show the result of the fitting
showFitting=1;

%
nDiffs=1;
immPop=0;
diffModel='unconfined';
pStart=[];

% filters to use when calculating the correlation functions
% 1:
% 2:
% 3:
% 4: non-overlapping time lags;
% 5:
sticsFilters=[5,2];

% which fitting parameters for the correlation gaussian fit?
% 1: amplitude
% 2: offset
% 3: x center
% 4: y center
% 5: x width
% 6: y width
% 7: angle
% 156 is independently variable amplitude and two widths. with offset 0 and
% specified rotation angle. 123456 is a translatable asymmetric gaussian
% with variable offset, amplitude and fixed (specified) angle.
gfittype=123456;
widths=find(num2str(gfittype)=='5'|num2str(gfittype)=='6');

% parameters for phase mask finding. dilate factor, low thresh, high thres,
% autofill, min area, max area
phaseParams=[2,0,.9,1,500,1000];

% maximum frame separation to consider
maxTau=4;

% magnification of the microscope in micrometers/pixel.
mpp=.049;

% time from beginning of one frame to the beginning of the next in seconds
intTime=.04;

tau=(1:maxTau)*intTime;

%% locate and prepare data
[datalist,dataloc,findex]=uigetfile([mainfold filesep ...
    '*.nd2;*.tif*;*.bin;*.xls*'],...
    'Matlab Files','multiselect','on');

if findex==0
    fprintf('no data selected\n')
    return
end

% make the output from uigetfile uniform regardless of number of files
% selected
if ~iscell(datalist); datalist={datalist}; end
for ii=1:numel(datalist); datalist{ii}=[dataloc datalist{ii}]; end
[dlocs,dnames,dexts]=cellfun(@fileparts,datalist,'uniformoutput',false);

% if a spreadsheet describing the data is selected, get the names from it
infomat=zeros(numel(dlocs),5);
if all(cellfun(@(x)~isempty(regexp(x,regexptranslate('wildcard','*xls*'),...
        'once')),dexts))
    [infomat,textdat]=xlsread(fullfile(dlocs{1},[dnames{1},dexts{ii}]));
    textdat=textdat(2:end,1:3);
    
    dlocs=textdat(:,2);
    dnames=textdat(:,1);
end

%% convert data to binary file if calculating correlation function and bin 
% doesn't exist
if yesStics==1
    for ii=1:numel(dnames)
        if strcmp(dexts,'.nd2')||strcmp(dexts,'.tif')||strcmp(dexts,'.tiff')
            writebin(fullfile(dlocs{ii},[dnames{ii},dexts]));
        end
    end
end

% Choose the functions for fitting
[corrFun,msdFun,pStart,corrLB,corrUB]=sticsFunFinder(nDiffs,immPop,diffModel,pStart);

% gooooooooooood luck figuring this out. linCell linearizes a cell array
% into a single 1D vector. linCell and longmsd both have to be aux
% functions in the code that fHandle and eHandle are used in.
fHandle=@(p,tau,corrDom,corrVal)nanLinCell(cellfun(@(x,y)x-y,...
    cellfun(@(x,y)corrFun(x,y,p),corrDom,num2cell(msdFun(tau,p),1),...
    'uniformoutput',0),corrVal,'uniformoutput',0));
eHandle=@(p,tau,corrDom)cellfun(@(x,y)corrFun(x,y,p),...
    corrDom,num2cell(msdFun(tau,p),1),'uniformoutput',0);

%% WRITE PHASEMASKS FILE FOR ALL MOVIES
if phasemasks&&yesStics&&checkvals
    if all(cellfun(@(x) regexp(x,regexptranslate('wildcard','*xls*')),dexts))
        plocs=dlocs;
        pnames=dnames;
        pext='.tif';
        
    else
        display('Select the phase images.')
        [phaselist,phaselistloc,findex]=uigetfile([mainfold filesep...
            '*.nd2;*.tif*;*.bin'],'Matlab Files','multiselect','on');
        if findex==0
            fprintf('no phase images selected. try again.\n')
            return
        end
        
        if ~iscell(phaselist); phaselist={phaselist}; end
        for ii=1:numel(phaselist); phaselist{ii}=[phaselistloc phaselist{ii}]; end
        [plocs,pnames,pext]=cellfun(@fileparts,phaselist,'uniformoutput',false);
    end
    
    for ii=1:numel(dnames)
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        if strcmp(pext,'.tif')||strcmp(pext,'.tiff')
            img=imread(fullfile(plocs{ii},[pnames{ii},pext]),1,...
                'Info',imfinfo(fullfile(plocs{ii},[pnames{ii},pext])));
            
        elseif strcmp(pext,'.nd2')
            vidid=bfGetReader(fullfile(plocs{ii},[pnames{ii},pext]));
            img=bfGetPlane(vidid,ii);
        end
        
        mnamelist=who(m); mnamelist=mnamelist(:);
        if any(cellfun(@strcmp,mnamelist,repmat({'phaseparams'},...
                [numel(mnamelist),1])))
            phaseParams=m.phaseparams;
        end
        
        fprintf(['phasing file named: ' dnames{ii} '.\n'])
        counter=0;
        while counter<2
            
            if exist('yn','var')&&~isempty(yn)
                phaseParams=yn;
                counter=0;
            else
                counter=counter+1;
            end
            
            phasemask=valley(img(:,:,1),phaseParams,checkvals);
            
            ppar.dilate_factor=phaseParams(1);
            ppar.low_threshold=phaseParams(2);
            ppar.high_threshold=phaseParams(3);
            ppar.autofill_bool=phaseParams(4);
            ppar.min_area=phaseParams(5);
            ppar.max_area=phaseParams(6);
            
            display(ppar)
            yn=input(['input new phase image parameters?\n enter for ''no'', '...
                'the parameters for ''yes''. \n two sequential ''nos'' means '...
                'it''s good to go.\n']);
        end
        if ~skipselect
            [~,goodcells]=select_cells(phasemask,img);
            phasemask(logical(-1*(ismember(phasemask,goodcells)-1)))=0;
        end
        m.phasemask=phasemask;
        m.phaseparams=phaseParams;
    end
elseif yesStics&&checkvals
    for ii=1:numel(dnames)
        [~,~,sz]=bingetframes([fullfile(dlocs{ii},dnames{ii}),'.bin'],1,[]);
        m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],'Writable',true);
        
        phasemask=ones(sz);
        m.phasemask=phasemask;
    end
end

if checkvals==1
    return          % no need to run any more code
end

if yesStics
    %% perform STICS
    for jj=1:numel(dnames)  % Loop each movie for sticsing
        
        % Display movie folder counter
        fprintf(['Working on this movie: ',dnames{jj},'\n'])
        
        m=matfile([fullfile(dlocs{jj},dnames{jj}),'_analysis.mat'],'Writable',true);
        phmask=m.phasemask;
        ncells=unique(phmask); ncells(ncells==0)=[];
        
        timecorrsave=cell(1,sum(ncells>0));
        thet=zeros(1,sum(ncells>0));
        leng=zeros(1,sum(ncells>0));
        nframes=zeros(1,sum(ncells>0));
        sphmask=cell(1,sum(ncells>0));
        counter=0;
        for ii=ncells(:)'
            counter=counter+1;
            
            sphmask{counter}=phmask==ii;
            
            roiinds=zeros(sum(double(sphmask{counter}(:))),2);
            [roiinds(:,1),roiinds(:,2)]=find(sphmask{counter});
            roi=[min(roiinds(:,1)),min(roiinds(:,2)),...
                max(roiinds(:,1)),max(roiinds(:,2))];
            
            % use phmask to assign the correct portion of phmask to another
            % variable, pm2
            sphmask{counter}=sphmask{counter}(roi(1):roi(3),roi(2):roi(4));
            
            [~,n]=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],1,roi);
            
            % max working data set size: 1 gigabyte/(8 bytes/double)
            maxdoubles=1e8/8;
            
            % movie is quadroupled in size due to fft padding
            potsize=n*(roi(3)-roi(1))*(roi(4)-roi(2))*4;
            
            % number of movie subdivisions
            nbins=ceil(potsize/maxdoubles);
            
            % number of frames in each subdivision
            binsize=floor(n/nbins);
            
            fprintf(['cross-correlating cell ', num2str(counter), ' of '...
                num2str(max(ncells)), '.\n'])
            if nbins>1
                
                timecorr=cell(1,nbins);
                for kk=1:nbins
                    fnums=1+binsize*(kk-1):binsize*kk;
                    v=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],...
                        fnums,roi);
                    timecorr{kk}=stics3(double(v),sphmask{counter},maxTau,...
                        sticsFilters);
                end
                
                timecorr=mean(cat(4,timecorr{:}),4);
            else
                
                v=bingetframes([fullfile(dlocs{jj},dnames{jj}),'.bin'],[],roi);
                timecorr=stics3(double(v),imdilate(sphmask{counter},ones(2)),...
                    maxTau,sticsFilters);
            end
            
            % use cyldist to find the cell's orientation and length
            [~,p,~,l]=cyldist(roiinds);
            thet(counter)=-atan(p(1)/p(2));
            leng(counter)=l-phaseParams(1)*2*mpp;
            
            timecorrsave{counter}=timecorr;
            nframes(counter)=n;
        end
        
        % WRITE RESULTS TO ANALYSIS FILE
        m.tcorr=timecorrsave;
        m.thet=thet;
        m.leng=leng;
        m.sticsfilters=sticsFilters;
        m.sphmask=sphmask;
        m.nframes=nframes;
    end
end
counter=0;
%% fit the correlations
if yesGauss
    for ii=1:numel(dnames)                % loop through movies
        if localFit
            m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],...
                'Writable',true);
            fprintf(['Fitting the correlations from this movie: ',dnames{ii},'\n'])
            
            tcorr=m.tcorr;
            thet=m.thet;
            sphmask=m.sphmask;
            sphmask=sphmask(cellfun(@(x)~isempty(x),sphmask));
            
            fCoeff=zeros(numel(tcorr),maxTau,numel(num2str(gfittype)));
            p95=fCoeff;
            nrmse=zeros(numel(tcorr),maxTau);
            for jj=find(cellfun(@(x)~isempty(x),tcorr)) % loop through cells in a movie
                counter=counter+1;
                
                % pad the phase mask to make it the size of the correlation func
                %         phmask=m.phasemask;
                bphmask=padarray(sphmask{jj},floor((size(tcorr{jj}(:,:,1))-...
                    size(sphmask{jj}))/2));
                if size(bphmask,1)<size(tcorr{jj},1)
                    bphmask=padarray(bphmask,[1,0],'pre');
                end
                if size(bphmask,2)<size(tcorr{jj},2)
                    bphmask=padarray(bphmask,[0,1],'pre');
                end
                
                res=zeros([size(tcorr{jj}(:,:,1)),maxTau]);
                for kk=1:size(tcorr{jj},3)-1  % loop through time lags
                    workingCorr=tcorr{jj}(:,:,kk+1);
                    workingCorr(~bphmask)=nan;
                    
                    [fCoeff(jj,kk,:),p95(jj,kk,:),nrmse(jj,kk),res(:,:,kk)]=...
                        gaussfit(workingCorr,gfittype,thet(jj));
                    
                end
                ampandwidths=find(num2str(gfittype)=='1'|num2str(gfittype)=='5'|...
                    num2str(gfittype)=='6');
                fCoeff(jj,:,ampandwidths)=abs(fCoeff(jj,:,ampandwidths));
                
                if showFitting
                    im1=tcorr{jj}(:,:,2); im1(~bphmask)=nan;
                    im2=tcorr{jj}(:,:,maxTau); im2(~bphmask)=nan;
                    r1=res(:,:,1); r1(~bphmask)=nan;
                    r2=res(:,:,maxTau); r2(~bphmask)=nan;
                    
                    subplot(321); pcolor(im1);
                    title('data')
                    axis image; shading flat
                    cl=get(gca,'clim');
                    
                    subplot(322); pcolor(r1);
                    title('residuals')
                    axis image; shading flat
                    cl=cl-min(cl)+nanmin(r1(:));
                    set(gca,'clim',cl);
                    
                    subplot(323); pcolor(im2);
                    title('data')
                    axis image; shading flat
                    cl=cl-min(cl)+nanmin(im2(:));
                    set(gca,'clim',cl);
                    
                    subplot(324); pcolor(r2);
                    title('residuals')
                    axis image; shading flat
                    cl=cl-min(cl)+nanmin(r2(:));
                    set(gca,'clim',cl);
                    
                    subplot(3,2,5:6)
                    scatter(1:size(fCoeff,2),...
                        squeeze(fCoeff(jj,:,6).^2*2*.049^2),'fill')
                    axis tight
                end
            end
            
        elseif ~localFit
            m=matfile([fullfile(dlocs{ii},dnames{ii}),'_analysis.mat'],...
                'Writable',true);
            fprintf(['Fitting the correlations from this movie: ',dnames{ii},'\n'])
            
            tcorr=m.tcorr;
            thet=m.thet;
            sphmask=m.sphmask;
            sphmask=sphmask(cellfun(@(x)~isempty(x),sphmask));
            
            nrmse=zeros(numel(tcorr),maxTau);
            for jj=find(cellfun(@(x)~isempty(x),tcorr)) % loop through cells in a movie
                counter=counter+1;
                
                % pad the phase mask to make it the size of the correlation func
                %         phmask=m.phasemask;
                bphmask=padarray(sphmask{jj},floor((size(tcorr{jj}(:,:,1))-...
                    size(sphmask{jj}))/2));
                if size(bphmask,1)<size(tcorr{jj},1)
                    bphmask=padarray(bphmask,[1,0],'pre');
                end
                if size(bphmask,2)<size(tcorr{jj},2)
                    bphmask=padarray(bphmask,[0,1],'pre');
                end
                
                workingCorr=tcorr{jj}(:,:,2:maxTau+1);
                workingCorr(~bphmask(:,:,ones(1,maxTau)))=nan;
                wSize=size(workingCorr);

                workingCorr=squeeze(mat2cell(workingCorr,wSize(1),wSize(2),...
                    ones(1,maxTau)))';
 
                [X,Y]=ndgrid(linspace(-wSize(1)/2,wSize(1)/2,wSize(1)),...
                    linspace(-wSize(1)/2,wSize(1)/2,wSize(2)));
                corrDomA=cat(3,X,Y);
                corrDomA(~bphmask(:,:,[1,1]))=nan;
                corrDom={corrDomA}; corrDom=corrDom(1,ones(1,maxTau));
                
                fCoeff(:,counter,jj)=lsqnonlin(@(p)fHandle(p,tau,corrDom,...
                    workingCorr),pStart,corrLB,corrUB);

                %% Plot the results
                predCorr=eHandle(fCoeff,tau,corrDom);
                startCorr=eHandle(pStart,tau,corrDom);
                residCorr=cellfun(@(x,y)x-y,workingCorr,predCorr,'uniformoutput',0);
                
                if showFitting
                    for kk=1:maxTau
                        subplot(3,maxTau,kk); pcolor(workingCorr{kk});
                        axis image; shading flat
                        cl=get(gca,'clim');
                        
                        subplot(3,maxTau,kk+maxTau); pcolor(residCorr{kk});
                        axis image; shading flat
                        cl=cl-min(cl)+nanmin(nanmin(residCorr{kk}));
                        set(gca,'clim',cl);
                        
                        subplot(3,maxTau,kk+2*maxTau); pcolor(startCorr{kk});
                        axis image; shading flat
                        cl=cl-min(cl)+nanmin(nanmin(startCorr{kk}));
                        set(gca,'clim',cl);
                    end
                end
            end
        end
        m.fCoeff=fCoeff*mpp.^2;
    end
end

if ~localFit
    output=fCoeff*mpp^2;
    return
end

%% select data

% load data
msds=cell(1,numel(dlocs));
cellLengths=cell(1,numel(dlocs));
infomatCell=cell(1,numel(dlocs));
nrmse=cell(1,numel(dlocs));
p95=cell(1,numel(dlocs));
nframes=cell(1,numel(dlocs));

fprintf('loading msd data\n')
parfor jj=1:numel(dlocs)
    m=matfile([fullfile(dlocs{jj},dnames{jj}),'_analysis.mat']);
    msds{jj}=m.pgauss.^2*2*mpp^2;
    cellLengths{jj}=m.leng;
    nrmse{jj}=m.nrmse;
    p95{jj}=m.p95;
    nframes{jj}=m.nframes;
    infomatCell{jj}=infomat(jj*ones(1,size(nframes{jj},2)),:);
    
end

cellLengths=cat(2,cellLengths{:})';
nrmse=cat(1,nrmse{:});
p95=cat(1,p95{:});
nframes=cat(2,nframes{:})';

infomat=cat(1,infomatCell{:});
infomat=cat(2,infomat,cellLengths);

msds=cat(1,msds{:});

%% fit msds
maxTau=15;
msdp=zeros(size(infomat,1),3,1);
dest=zeros(size(infomat,1),1);
msdp95=zeros(size(infomat,1),3,1);
fprintf('fitting msds\n')
for jj=1:floor(size(infomat,1))
    intTime=infomat(jj,1)*.001;
    tau=(1:maxTau)*intTime;
    pstart=[cellLengths(jj)-.5,5,.0384];
    Y=msds(jj,1:maxTau,6);
    
    lb=[0,0,0]; ub=[100,200,100];
    
    %     mdl=fitnlm(tau,y,@confmodel,pstart);
    %     msdp(jj,:,1)=abs(mdl.Coefficients{:,1});
    
    [X,~,resid,~,~,~,jaco]=lsqcurvefit(@longmsd1d,pstart,tau,Y,lb,ub);
    
    msdp95(jj,:)=diff(nlparci(X,resid,'jacobian',jaco),1,2);
    
    msdp(jj,:)=X;
    
    dest(jj,1)=(Y(2)-Y(1))/4/intTime;
    %     msdp95(jj,:,1)=diff(coefCI(mdl),[],2);
    
end

output.msdp=msdp;
output.dEstimate=dest;
output.msdp95=msdp95;
output.nframes=nframes;
output.cellLengths=cellLengths;
output.inttime=infomat(:,1)*.001;
output.infomat=infomat;
end

% square confinement model
function z=longmsd1d(p,x)
% global camerasd
d=p(1); l=p(2); tau=x;

summedterm=@(t,d,l,n)1/n^4*exp(-(n*pi/l).^2*d*t);

temp=eps*ones(size(tau));
for ii=1:2:2*400-1
    s=summedterm(tau,d,l,ii);
    if sum(s./temp)<1e-10
        break
    end
    temp=temp+s;
end
z=l^2/6*(1-96/pi^4*temp)+p(3);
end

function [cell_xy,good_cell]=select_cells(PhaseMask,im)

% Let user click on the phase mask of cell images to decide which cell to
% analyze subsequently

% INPUTS:

% PhaseMask: Numbered phase mask of cell images returned by the 'valley.m'
% code or any other cell segmentation code

% OUTPUTS:

% cell_xy: The x- and y-coordinates clicked on by the users to select cells

% good_cell: The indices of chosen cells

cell_fig_h=figure;
subplot(121)
imshow(im,[])
subplot(122)
imshow(PhaseMask~=0,[],'Border','Tight');
title(sprintf(['Left-click on cells to be analyzed.\n Press ''Enter'' to proceed.\n'...
    ' Or click return without clicking\n on any cell to analyze ALL of them.']))

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
    [curr_cell_x,curr_cell_y,key_input]=jsbginput(1);
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

function out=nanLinCell(in)
out=cat(1,in{:});
out=out(~isnan(out));
end