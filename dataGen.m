function [v, simProps] = dataGen(varargin)
%
% NAME:
%       dataGen
% PURPOSE:
%       Generates space and time resolved single-molecule imaging data of a
% 			single diffuser.
% CATEGORY:
%       Data Simulation
% CALLING SEQUENCE:
%       [v, simProps] = dataGen(boundaryCondition, blurFlag);
% INPUTS:
%       varargin:  use paired inputs to set the property (input 1) to the
%           value (input 2) desired.
%
%       Properties:             Descriptions:
%
%       D                       diffusion coefficient in microns^2/s
%
%       tFrame                  image frame integration time in seconds
%
%       pixSize                 Pixel size in micrometers
%
%       psfSize                 standard deviation (width) of the microscope's
%                               point spread function in micrometers. i.e.
%                               FWHM = sqrt(2*log(2)) * psfSize
%
%       celSize                 1x2 vector: [width (diameter), height
%                               (length)] of the bounding box (cylinder)
%                               depending on the confinement condition
%
%       nFrames                 number of frames to be simulated
%
%       SNR                     signal to noise ratio (ratio of maximum
%                               signal amplitude to standard deviation of
%                               background noise)
%
%       confBool                1 for confined to the interior of a
%                               cylinder, and 0 for free diffusion.
%       
%       blurFlag                1 for blurry motion, 0 for 'stroboscopic
%                               illumination'
%
% OUTPUTS:
%       v:            		simulated image time sequence
%       simProps:        	properties of the simulation in the Matlab structure format
% PROCEDURE:
%       1. Simulate molecular trajectory
%       2. Evaluate pixel intensities
%		3. Add noise
% MODIFICATION HISTORY:
%       Written by David J. Rowland, The University of Michigan, 3/16.
% NOTES:
%       This code 'dataGen.m' should be considered 'freeware'- and may be
%       distributed freely in its original form when properly attributed.
%
%       For testing purposes, this line makes a movie of confined, blurry
%       diffusion.
%
%       v=dataGen('confined', 1)


%% Default simulation parameters
D = .01;                % diffusion coefficient in microns^2/s
tFrame = .05;			% frame integration time in seconds
pixSize = .049;			% width of pixels in microns
psfSize = .098;			% s.d. of the psf in microns
celSize = [1,3];        % [width, length] of confinement cylinder in microns
nFrames = 1e3;			% number of frames in the simulated movie
SNR = 20;				% signal to noise ratio for added white noise (0:inf)
confBool = 1;           % 'confined' or 'unconfined'
blurFlag = 1;           % include blur subframes or not

% algorithmically-determined image size designed to disallow edge effects
imSize = ceil(celSize/pixSize+4*ceil(psfSize/pixSize));

% check variable size
if prod([imSize(1:2),nFrames])*8/1e9>1
    warning('video is over a GB. manually pass this block if you wish to continue')
end

% minimum increment of speed in units of microns^2/subframe
dtRef = 0.0001;

% initialize simulation properties structure
simProps.D = D;
simProps.tFrame = tFrame;
simProps.pixSize = pixSize;
simProps.psfSize = psfSize;
simProps.celSize = celSize;
simProps.nFrames = nFrames;
simProps.SNR = SNR;
simProps.blurFlag = blurFlag;
simProps.nSubs = [];
simProps.dtRef = dtRef;
simProps.confBool = confBool;

% if any sim parameters are included as inputs, change the simulation
% parameters mentioned
if ~rem(nargin,2)
    fNames=fieldnames(simProps);
    for ii=1:2:nargin
        whichField = strcmp(fNames,varargin{ii});
        
        if all(~whichField)
            warning('Check spelling. Parameter change may have not occurred')
        end
        
        eval([fNames{whichField} ' = varargin{ii+1}'])
        eval(['simProps.' fNames{whichField} ' = ' fNames{whichField},';'])
    end
    
elseif ~rem(nargin,1)
    warning('use paired inputs')
    v=[];
    return
end

% number of subframes required
nSubs = ceil(D*tFrame/dtRef);

% update the value of the diffusion coefficient since rounding may change it.
D = nSubs*dtRef/tFrame;

%% Trajectory generation
if confBool
    % confined particle trajectory
    mLocs = zeros(3, nFrames);
    for ii = 1 : nFrames*nSubs-1    % this loop can be compiled to mex64 to increase its speed
        % three 1d steps pulled from normal distribution with variance 2*dtRef
        step = sqrt(2*dtRef) * randn(3,1);
        
        candPos = mLocs(:,ii) + step;
        prevPos = mLocs(:,ii);
        
        % the ordering in the celSize vector matters because of this line:
        r = celSize(1)/2;
        
        % if the candidate position is outside of the cylinder in the x/z dimensions, reflect the step
        % against the inside of the cylinder. path length is preserved.
        if sqrt(sum(candPos([1,3]).^2)) > celSize(1)/2
            m = (candPos(3)-prevPos(3)) / (candPos(1)-prevPos(1));
            b = candPos(3)-m*candPos(1);
            
            xi=[(-m*b+sqrt(-b^2+r^2+m^2*r^2))/(1+m^2),...
                (-m*b-sqrt(-b^2+r^2+m^2*r^2))/(1+m^2)];
            yi=m*xi+b;
            
            % there are two solutions. the one closest to the candidate position is chosen. the farther
            % one is on the other side of the cell.
            whichone = (xi-candPos(1)).^2 + (yi-candPos(3)).^2 + (xi-prevPos(1)).^2 + (yi-prevPos(3)).^2;
            xip = [xi(find(whichone == min(whichone))),yi(find(whichone == min(whichone)))];
            
            normv = -xip/sqrt(sum(xip.^2));
            l = sqrt(sum((candPos([1,3])'-xip).^2));
            pf = 2*sum((prevPos([1,3])'-xip).*normv)*normv-(prevPos([1,3])'-xip);
            pf = pf/sqrt(sum(pf.^2))*l+xip;
            
            % replace x/z components of the position with the reflected x/z components
            out=candPos;
            out([1,3])=pf;
            
            candPos=out;
        end
        
        % if the candidate position is outside of the cylinder in the y dimension (cell's long axis)
        if candPos(2) < -celSize(2)/2
            candPos(2) = 2*-celSize(2)/2 - candPos(2);
        end
        if candPos(2) > celSize(2)/2
            candPos(2) = 2*celSize(2)/2 - candPos(2);
        end
        
        mLocs(:, ii+1) = candPos;
    end
else
    
    % unconfined particle trajectory
    mLocs=cumsum(sqrt(2*dtRef) * randn(3,nFrames*nSubs),2);
end

%% movie generation

if blurFlag
    % arrange tracks for subframe averaging
    tr_x = zeros(nSubs, nFrames);
    tr_y = zeros(nSubs, nFrames);
    for ii = 1:nFrames
        tr_x(:, ii) = mLocs(1, 1+(ii-1)*nSubs : ii*nSubs);
        tr_y(:, ii) = mLocs(2, 1+(ii-1)*nSubs : ii*nSubs);
    end
else
    % just use the first subframe from each frame
    tr_x = mLocs(1, 1:nSubs:end);
    tr_y = mLocs(2, 1:nSubs:end);
end

% shift to positive values
tr_x=tr_x+celSize(1)/2;
tr_y=tr_y+celSize(2)/2;

% noiseless pixel intensities
v=zeros([imSize(1:2),nFrames]); cx = 0; cy = 0;
for ii=-2*psfSize:pixSize:celSize(1)+2*psfSize        % x pixel locations with padding
    cx = cx+1;
    for jj=-2*psfSize:pixSize:celSize(2)+2*psfSize    % y pixel locations with padding
        cy = cy+1;
        
        % symmeteric gaussian function approximation of Airy Disk
        v(cx,cy,:) = mean(exp(-((ii-tr_x).^2+(jj-tr_y).^2)/2/psfSize^2),1);
    end
    cy = 0;
end

% add white noise
v = v + 1/SNR*randn(size(v));

simProps.D = D;
simProps.tFrame = tFrame;
simProps.pixSize = pixSize;
simProps.psfSize = psfSize;
simProps.celSize = celSize;
simProps.nFrames = nFrames;
simProps.SNR = SNR;
simProps.blurFlag = blurFlag;
simProps.nSubs = nSubs;
simProps.dtRef = dtRef;
simProps.confBool = confBool;
end