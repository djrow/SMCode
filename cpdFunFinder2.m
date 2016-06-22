function outStruct = cpdFunFinder2(anProp)

dim = anProp.dim;
nDiffs = anProp.nMobile;
immBool = anProp.immBool;
confBool = anProp.confBool;
globBool = anProp.globBool;

% partial 2d cpd function
c2=@(x,y,p)p*exp(-x./y);

% partial 1d cpd function
c1=@(x,y,p)p*erf(sqrt(x./2/y));

% 2d confined msd function
m2=@(t,p)4*p(1)*t+p(2);

% 1d unconfined msd function
m1=@(t,p)2*p(1)*t+p(2);

pStart = [.1, 0, .05,eps, .01,eps, .0025,eps, ...
    .25, .25, .25, 1];

LB=zeros(1,numel(pStart));
UB=inf(1,numel(pStart));
UB(9:11) = 1;

if dim == 2 && immBool == 0 && globBool == 1
    if ~confBool
        switch nDiffs
            case 1
                msdFun=@(tau,p) ...
                    m2(tau,p([1,2]));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),1);
                pID = 1:2;
                
            case 2
                msdFun=@(tau,p)cat(2,...
                    m2(tau,p([1,2])),...
                    m2(tau,p([3,2])));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),p(4))-...
                    c2(x,y(2),1-p(4));
                pID = [1:3,9];
                
            case 3
                msdFun=@(tau,p)cat(2,...
                    m2(tau,p([1,2])),...
                    m2(tau,p([3,2])),...
                    m2(tau,p([4,2])));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),p(5))-...
                    c2(x,y(2),p(6))-...
                    c2(x,y(3),1-p(5)-p(6));
                pID = [1:3,5,9:10];
                
            case 4
                msdFun=@(tau,p)cat(2,...
                    m2(tau,p([1,2])),...
                    m2(tau,p([3,2])),...
                    m2(tau,p([4,2])),...
                    m2(tau,p([5,2])));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),p(6))-...
                    c2(x,y(2),p(7))-...
                    c2(x,y(3),p(8))-...
                    c2(x,y(4),1-p(6)-p(7)-p(8));
                pID = [1:3,5,7,9:11];
        end
    else
        switch nDiffs
            case 1
                msdFun=@(tau,p) ...
                    confMSD2(tau,p([1,2,3]));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),1);
                pID = [1:2,12];
                
            case 2
                msdFun=@(tau,p)cat(2,...
                    confMSD2(tau,p([1,2,5])),...
                    confMSD2(tau,p([3,2,5])));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),p(4))-...
                    c2(x,y(2),1-p(4));
                pID = [1:3,9,12];
                
            case 3
                msdFun=@(tau,p)cat(2,...
                    confMSD2(tau,p([1,2,7])),...
                    confMSD2(tau,p([3,2,7])),...
                    confMSD2(tau,p([4,2,7])));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),p(5))-...
                    c2(x,y(2),p(6))-...
                    c2(x,y(3),1-p(5)-p(6));
                pID = [1:3,5,9:10,12];
                
            case 4
                msdFun=@(tau,p)cat(2,...
                    confMSD2(tau,p([1,2,9])),...
                    confMSD2(tau,p([3,2,9])),...
                    confMSD2(tau,p([4,2,9])),...
                    confMSD2(tau,p([5,2,9])));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),p(6))-...
                    c2(x,y(2),p(7))-...
                    c2(x,y(3),p(8))-...
                    c2(x,y(4),1-p(6)-p(7)-p(8));
                pID = [1:3,5,7,9:12];
        end
    end
    pStart={pStart(pID)};
    bounds=[{LB(pID)},{UB(pID)}];
    dID = find(ismember(pID,1:2:7));
    aID = find(ismember(pID,9:11));
    
    
elseif ~globBool && ~confBool
    msdFun = @(tau,p)m2(tau,p);
    msdStart = pStart([1,2]);
    msdLB = [0,-inf];
    msdUB = [inf,inf];
    switch nDiffs
        case 1
            cpdFun = @(x,p)1-...
                c2(x,p(1),1);
            cpdStart = pStart(1);
            cpdLB = 0;
            cpdUB = inf;
        case 2
            cpdFun = @(x,p)1-...
                c2(x,p(1),p(3))-...
                c2(x,p(2),1-p(3));
            cpdStart = pStart([1,3,9]);
            cpdLB = [0,0,0];
            cpdUB = [inf,inf,1];
    end
    pStart = [{cpdStart},{msdStart}];
    bounds = [{cpdLB},{cpdUB},{msdLB},{msdUB}];
    dID = nan;
    aID = nan;
else
    warning('unsupported parameters')
    outStruct = [];
end
outStruct.cpdFun = cpdFun;
outStruct.msdFun = msdFun;
outStruct.pStart = pStart;
outStruct.bounds = bounds;
outStruct.dID = dID;
outStruct.aID = aID;
end

%% 2d square confinement model
function z=confMSD2(tau,p)
% global camerasd
d=p(1); l=p(2);

summedterm=@(t,d,l,n)1/n^4*exp(-(n*pi/l).^2*d*t);

temp=eps*ones(size(tau));
for ii=1:2:2*400-1
    s=summedterm(tau,d,l,ii);
    if sum(s./temp)<1e-10
        break
    end
    temp=temp+s;
end
z=l^2/3*(1-96/pi^4*temp)+p(3);
end