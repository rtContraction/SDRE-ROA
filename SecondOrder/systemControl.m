function [u,E,dE]=systemControl(x,xr,ur,sdre,sys)
% ur,xr= 0 for the equilibrium point
option=sdre.comp;
    switch option
        case 1  % CLF
            [u,E,dE]=u_CLF(x);  
        case 2  % SDRE
            [u,E,dE]=u_SDRE(x,sys,sdre);
    end
end

function [u,E,dE]=u_CLF(x)
x1=x(1,:);
x2=x(2,:);
k1=1;
k2=1;
z=x2+(x1^2+k1*x1);
u=2*k1*x1^2-((k1/k2)-k1^2)*x1-(2*x1+k1+1)*z;
E=.5*k1*x1^2+.5*k2*z^2;
dE=(k1*x1)*(-k1*x1+z)+k2*z*(-2*k1*x1^2-k1^2*x1+2*x1*z+k1*z+u);
end

function [u,E,dE]=u_SDRE(x,sys,sdre)
A=sys.A(x);            % A(x)*x=f(x)
B=sys.B;
[K,~,~] = lqr(A,B,sdre.Q,sdre.R); 
u=-K*x;                     % control law
f_x=(A-B*K);                % closed-loop system
M=sdre.M(x);                % contraction Metric
Md=sdre.Md([x;f_x(:)]);     % dM(x)  
E=x'*M*x;                   % distance Function
Jf =[2*x(1)  1;0 0];
F=Jf-B*K;                   % differential system
dE=x'*(Md+M*F+F.'*M)*x;     % time derivative of square distance
end
