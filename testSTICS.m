function [v,caseInfo]=testSTICS(dType,snr)
% this code generates movies that can be used to test the STICS code
% it first uses testCPD to generate trajectories.

% types of diffusion simulatable
diffusionTitles{1}='one mobile';
diffusionTitles{3}='one mobile, one immobile';
diffusionTitles{5}='two mobile';
diffusionTitles{7}='two mobile, one immobile';
diffusionTitles{9}='one confined mobile';
diffusionTitles{11}='one confined mobile, one immobile';
diffusionTitles{13}='two confined mobile';
diffusionTitles{15}='two confined mobile, one immobile';

%% testCPD trajectories are in pixels
[tr,caseInfo]=testCPD(dType);

mpp=caseInfo.mag;

% parameters specific to imaging the trajectories
sdCamera=2*caseInfo.mag;
caseInfo.sdCamera=sdCamera;
% snr=inf;
caseInfo.imagingSNR=snr;

xqPixels=(((min(tr(:,4))-sdCamera*4)/mpp):1:...
    (max(tr(:,4))+sdCamera*4)/mpp)*mpp;
yqPixels=(((min(tr(:,5))-sdCamera*4)/mpp):1:...
    (max(tr(:,5))+sdCamera*4)/mpp)*mpp;

if numel(xqPixels)>500
    warning('too large of a movie')
    return
end

%% find the total brightness of every pixel.
v=zeros(numel(xqPixels),numel(yqPixels),size(tr,1));
for ii=1:numel(xqPixels)
    for jj=1:numel(yqPixels)
        v(ii,jj,:)=exp(-((xqPixels(ii)-tr(:,4)).^2+...
            (yqPixels(jj)-tr(:,5)).^2)/2/sdCamera^2);
    end
end

%% noise
if ~isinf(snr)
    nterm=randn(size(v))/snr;
    
    nterm=nterm-min(nterm(:));
    v=v+nterm;
    nval=max(v(:));
    v=v/nval;
end
end