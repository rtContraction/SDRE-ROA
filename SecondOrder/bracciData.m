function [PclA,lyapFun] = bracciData(xrange,lineWidthColor,varargin)

%==================   System Properties  ==================================
b1=0;
b2=1;
B=[b1;b2];
n=2;
A_x = @(x) [x(1) 1; 0 0];
%==================   SDRE Parameters =====================================
R=1;
Q = diag([100,100]);
% ===========================  SDP Variables  =============================
A=A_x(zeros(2,1));
[~,gainSol,~]=icare(A,B,Q,R);  %Q(x3) ?, R(x3)?
A0=A - B*gainSol;
P0=lyap(A0',Q);
% ===========================  SDP Variables  =============================
z=sdpvar(n,1); % free variable for optimasation
lyapFun=matlabFunction(SDRE_SDP2Sym(z,z.'*P0*z),'Vars',...
{SDRE_SDP2Sym(z,z(1)),SDRE_SDP2Sym(z,z(2))});
PclA=[];
%==================  Gridding Terms =======================================
if nargin>2 % for analysis
nPoints=varargin{1};
xPoints1=nPoints(1)+1; 
xPoints2=nPoints(2)+1;
xP1=linspace(-xrange(1),xrange(1),xPoints1);
xP2=linspace(-xrange(2),xrange(2),xPoints2);
Acl=cell(xPoints1,xPoints2);
Fcl=cell(xPoints1,xPoints2);
Pcl=cell(xPoints1,xPoints2);
L=xPoints1*xPoints2; % Total number of points

% ===================  1. Gridding Points  ================================
barGUI = waitbar(1/L,'1','Name','SDRE Simulation running');               
for j=1:xPoints2
    for i=1:xPoints1
        x=[xP1(i);xP2(j)];
        %  ==================== SDRE ======================================
        [Psol,gainSol,eigSol]=icare(A_x(x),B,Q,R);  %Q(x3) ?, R(x3)?
        Acli=A_x(x) - B*gainSol;
        Acl{i,j}=Acli;
        % ===========  Lyapunov Condition  =============================
        F=2*x'*P0*Acli*x; % closed-loop system
        Fcl{i,j}=F;
        if F<0
            Pcl{i,j}=x;
        else
            Pcl{i,j}=[0;0];
        end
        waitbar((i+xPoints2*(j-1))/L,barGUI,...
            sprintf('Pointwise evaluation at:%0.5f,%0.5f',x(1),x(2)))
    end
end
waitbar(1,barGUI,'Checking Pointwise Conditions...');
close(barGUI)
cShade=[255,0,255]/255;
PclA=cell2mat(Pcl);
PPcl=reshape(PclA,[2 size(PclA,1)/2 * size(PclA,2)]);
plot(PPcl(1,:),PPcl(2,:),'Marker','x','LineStyle','none','Color',cShade)
end

[X,Y] = meshgrid(-xrange(1):0.05:xrange(1),-xrange(2):0.05:xrange(2));
Z=lyapFun(X,Y);
hold on
E0=65;
contour(X,Y,Z,[E0 E0],':m','linewidth',lineWidthColor);
end
