function [sticsFun,msdFun,pstart,cpdLB,cpdUB]=sticsFunFinder(nDiffs,immPop,dType,pstart)

switch immPop
    case 0 % no immobile term
        switch nDiffs
            case 1 % 1 diffuser
                switch dType
                    case 'unconfined' % 1 unconfined diffuser in 2D
                        sticsFun=@(x,y,p)p(3)+p(4)*exp((-(x(:,:,1)-p(5)).^2-...
                            (x(:,:,2)-p(6)).^2)/(2*y))./(2*pi*y);
                        msdFun=@(tau,p)2*p(1)*tau+p(2);
                        
                        if isempty(pstart)
                            pstart=[.1,.1,0,1,0,0];
                        end
                        cpdLB=[0,-inf,-inf,0,-inf,-inf];
                        cpdUB=inf(1,6);
                    case 'confined' % 1 confined diffuser in 2D
                        sticsFun=@(x,y,p)p(4)+p(5)*exp((-(x(:,:,1)-p(6)).^2-...
                            (x(:,:,2)-p(7)).^2)/(2*y))./(2*pi*y);
                        msdFun=@(tau,p)abs(longmsd1d(p(1:3),tau));
                        
                        if isempty(pstart)
                            pstart=[.1,1,.1,0,1,0,0];
                        end
                        cpdLB=[0,0,-inf,-inf,0,-inf,-inf];
                        cpdUB=inf(1,7);
                end
            case 2 % 2 diffusers
                switch dType
                    case 'unconfined' % 2 unconfined diffusers in 2D
                        sticsFun=@(x,y,p)p(5)+p(6)*...
                            exp((-(x(:,:,1)-p(7)).^2-(x(:,:,2)-p(8)).^2)/(2*y(1)))./(2*pi*y(1))+...
                            p(9)*...
                            exp((-(x(:,:,1)-p(7)).^2-(x(:,:,2)-p(8)).^2)/(2*y(2)))./(2*pi*y(2));
                        msdFun=@(tau,p)cat(1,2*p(1)*tau+p(2),2*p(3)*tau+p(4));
                        
                        if isempty(pstart)
                            pstart=[.1,.01,.001,.01,0,1,0,0,1];
                        end
                        cpdLB=[0,-inf,0,-inf,-inf,0,-inf,-inf,0];
                        cpdUB=inf(1,9);
                    case 'confined' % 2 confined diffusers in 2D
                        sticsFun=@(x,y,p)p(7)+p(8)*...
                            exp((-(x(:,:,1)-p(9)).^2-(x(:,:,2)-p(10)).^2)/(2*y(1)))./(2*pi*y(1))+...
                            p(11)*...
                            exp((-(x(:,:,1)-p(9)).^2-(x(:,:,2)-p(10)).^2)/(2*y(2)))./(2*pi*y(2));
                        msdFun=@(tau,p)cat(1,longmsd1d(p(1:3),tau),longmsd1d(p(4:6),tau));
                        
                        if isempty(pstart)
                            pstart=[.3,1,.01,.01,1,.01,0,1,0,0,1];
                        end
                        cpdLB=[0,0,-inf,0,0,-inf,-inf,0,-inf,-inf,0];
                        cpdUB=inf(1,11);
                end
        end
    case 1 % 1 immobile term
        switch nDiffs
            case 1 % 1 diffuser
                switch dType
                    case 'unconfined' % 1 unconfined diffuser with 1 immobile population in 2D
                        sticsFun=@(x,y,p)p(3)+p(4)*...
                            exp((-(x(:,:,1)-p(5)).^2-(x(:,:,2)-p(6)).^2)/(2*y(1)))./(2*pi*y(1))+...
                            p(7)*...
                            exp((-(x(:,:,1)-p(5)).^2-(x(:,:,2)-p(6)).^2)/(2*p(2)))./(2*pi*p(2));
                        msdFun=@(tau,p)2*p(1)*tau+p(2);
                        
                        if isempty(pstart)
                            pstart=[.1,.01,0,1,0,0,1];
                        end
                        cpdLB=[0,-inf,-inf,0,-inf,-inf,0];
                        cpdUB=inf(1,7);
                    case 'confined' % 1 confined diffuser with 1 immobile population in 2D
                        sticsFun=@(x,y,p)p(4)+p(5)*...
                            exp((-(x(:,:,1)-p(6)).^2-(x(:,:,2)-p(7)).^2)/(2*y(1)))./(2*pi*y(1))+...
                            p(8)*...
                            exp((-(x(:,:,1)-p(6)).^2-(x(:,:,2)-p(7)).^2)/(2*p(3)))./(2*pi*p(3));
                        msdFun=@(tau,p)longmsd1d(p(1:3),tau);
                        
                        if isempty(pstart)
                            pstart=[.1,1,.01,0,1,0,0,1];
                        end
                        cpdLB=[0,0,-inf,-inf,0,-inf,-inf,0];
                        cpdUB=inf(1,8);
                end
            case 2 % 2 diffusers
%                 switch dType
%                     case 'unconfined' % 2 unconfined diffusers with 1 immobile population in 2D
%                         sticsFun=@(x,y,p)1-p(5)*exp(-x./(y+p(7)))-...
%                             p(6)*exp(-x./(y+p(7)))-(1-p(5)-p(6))*exp(-x./p(7));
%                         msdFun=@(tau,p)cat(1,4*p(1)*tau+p(2),4*p(3)*tau+p(4));
%                         
%                         if isempty(pstart)
%                             pstart=[.1,eps,.01,eps,.3,.3,1e-3];
%                         end
%                         cpdLB=[0,-inf,0,-inf,0,0,0];
%                         cpdUB=[inf,inf,inf,inf,1,1,inf];
%                     case 'confined' % 2 confined diffusers with 1 immobile population in 2D
%                         
%                         sticsFun=@(x,y,p)1-p(7)*exp(-x./(y+p(9)))-...
%                             p(8)*exp(-x./(y+p(9)))-(1-p(7)-p(8))*exp(-x./p(9));
%                         msdFun=@(tau,p)cat(1,longmsd2d(p(1:3),tau),longmsd2d(p(4:6),tau));
%                         
%                         if isempty(pstart)
%                             pstart=[.1,.5,eps,.01,.5,eps,.3,.3,1e-3];
%                         end
%                         cpdLB=[0,0,-inf,0,0,-inf,0,0,0];
%                         cpdUB=[inf,inf,inf,inf,inf,inf,1,1,inf];
%                 end
        end
end


% % gooooooooooood luck figuring this out. linCell linearizes a cell array
% % into a single 1D vector. linCell and longmsd both have to be aux
% % functions in the code that fHandle and eHandle are used in.
% fHandle=@(p,tau,sqSteps,ranks)linCell(cellfun(@(x,y)x-y,...
%     cellfun(@(x,y)cpdFun(x,y,p),sqSteps,num2cell(msdFun(tau,p),1)',...
%     'uniformoutput',0),ranks,'uniformoutput',0));
% eHandle=@(p,tau,sqSteps)cellfun(@(x,y)cpdFun(x,y,p),...
%     sqSteps,num2cell(msdFun(tau,p))','uniformoutput',0);
end

%% square confinement model
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