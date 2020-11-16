function [u,E,dE]=systemControl(x,xr,ur,sdre,sys)
    option=sdre.comp;
    switch option
         case 1  % SDRE
            [u,E,dE]=u_LQR(x,xr,sys,sdre);   
        case 2  % SDRE
            [u,E,dE]=u_LQR(x,xr,sys,sdre);
        case 3  % SDRE
            [u,E,dE]=u_SDRE(x,xr,sys,sdre); 
    end

end

function [u,E,dE]=u_LQR(x,xr,sys,sdre)
[A, B] = systemSDC(x*0);
Q=diag([150,3,51,2,150,1]); % Q at x=0;
if sdre.comp==2
    gain=1e5;           % restricted controller
else
    gain=1e5-9.95e4;    % aggresive controller
end
R = eye(2).*gain; 
[P, K, ~] = icare(A,B,Q,R,[],[],[]); 
[A, B] = systemSDC(x); % real system
u=-K*x;
E=x'*P*x;
f_x=(A-B*K);
dE=2*x'*P*f_x*x;
end

function [u,E,dE]=u_SDRE(x,xr,sys,sdre)
[A, B] = systemSDC(x);
[Q, R] = systemWeightingMatrices(x);
[~, K, ~] = icare(A,B,Q,R,[],[],[]); 
u=-K*x;
M=sdre.M(x);
E=sdre.E(x(1),x(2),x(3),x(4),x(5),x(6));
A_cl=(A-B*K);
dA_dx_x = systemJacobianOpenLoop(x);
dBKx_dx = systemJacobianStateFeedback(x, K);                            
F = A + dA_dx_x - B*K - dBKx_dx;
Md=sdre.Md([x;A_cl(:)]);
dE=x'*(Md+M*F+F'*M)*x;

end
