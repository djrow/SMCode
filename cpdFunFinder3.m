function outStruct = cpdFunFinder3(anProp)
% global and local fitting functions for 2d unconfined diffusion

nDiffs = anProp.nMobile;
globBool = anProp.globBool;

% partial 2d cpd function
c2=@(x,y,p)p*exp(-x./y);

% 2d confined msd function
m2=@(t,p)4*p(1)*t+p(2);

pStart = [.1,.1,.01, .001, .0001, .00001, .2, .2, .2, .2, 1];

LB=zeros(1,numel(pStart));
UB=inf(1,numel(pStart));
UB(7:10) = 1;

switch globBool
    case 1      % global fitting (one fit for both CPD and MSD curves)
        switch nDiffs
            case 1      % one diffusion population
                msdFun=@(tau,p) ...
                    m2(tau,p([1,2]));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),1);
                pID = 1:2;
                
            case 2      % two diffusion populations
                msdFun=@(tau,p)cat(2,...
                    m2(tau,p([1,2])),...
                    m2(tau,p([3,2])));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),p(4))-...
                    c2(x,y(2),1-p(4));
                pID = [1:3,7];
                
            case 3      % three diffusion populations
                msdFun=@(tau,p)cat(2,...
                    m2(tau,p([1,2])),...
                    m2(tau,p([3,2])),...
                    m2(tau,p([4,2])));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),p(5))-...
                    c2(x,y(2),p(6))-...
                    c2(x,y(3),1-p(5)-p(6));
                pID = [1:4,7:8];
                
            case 4      % four diffusion populations
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
                pID = [1:5,7:9];
                
            case 5      % five diffusion populations
                msdFun=@(tau,p)cat(2,...
                    m2(tau,p([1,2])),...
                    m2(tau,p([3,2])),...
                    m2(tau,p([4,2])),...
                    m2(tau,p([5,2])),...
                    m2(tau,p([6,2])));
                cpdFun=@(x,y,p)1-...
                    c2(x,y(1),p(7))-...
                    c2(x,y(2),p(8))-...
                    c2(x,y(3),p(9))-...
                    c2(x,y(4),p(10))-...
                    c2(x,y(5),1-p(7)-p(8)-p(9)-p(10));
                pID = [1:10];
        end
        pStart={pStart(pID)};
        bounds=[{LB(pID)},{UB(pID)}];
        dID = find(ismember(pID,[1,3:6]));
        aID = find(ismember(pID,7:10));
        
    case 0      % local fitting (separate CPD and MSD fits)
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
                cpdStart = pStart([1,3,7]);
                cpdLB = [0,0,0];
                cpdUB = [inf,inf,1];
            case 3
                cpdFun = @(x,p)1-...
                    c2(x,p(1),p(4))-...
                    c2(x,p(2),p(5))-...
                    c2(x,p(3),1-p(4)-p(5));
                cpdStart = pStart([1,3,4,7,8]);
                cpdLB = zeros(1,5);
                cpdUB = [inf(1,3),1,1];
            case 4
                cpdFun = @(x,p)1-...
                    c2(x,p(1),p(5))-...
                    c2(x,p(2),p(6))-...
                    c2(x,p(3),p(7))-...
                    c2(x,p(4),1-p(5)-p(6)-p(7));
                cpdStart = pStart([1,3,4,5,7:9]);
                cpdLB = zeros(1,7);
                cpdUB = [inf(1,4),1,1,1];
        end
        pStart = [{cpdStart},{msdStart}];
        bounds = [{cpdLB},{cpdUB},{msdLB},{msdUB}];
        dID = nan;
        aID = nan;
end
outStruct.cpdFun = cpdFun;
outStruct.msdFun = msdFun;
outStruct.pStart = pStart;
outStruct.bounds = bounds;
outStruct.dID = dID;
outStruct.aID = aID;
end