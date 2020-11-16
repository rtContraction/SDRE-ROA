function Wd=SDRE_Md(n,x,f,W,varargin)
Wd=sdpvar(n);
dW=jacobian(W,x); %each layer is dW/dxj= dW11,1 dW12,1
if  nargin < 5
    % Wp= dW/dx *xp
    dW=shiftdim(dW,1); %Jacobian are rows dW11,1 dW11,2 ...
    for i=1:n
        for j=1:n
        Wd(i,j)=reshape(dW(j,:,i),[1 n])*f;
        end
    end
else
    % d(W(x)*x)
    for i=1:n
        for j=1:n
            Wd(i,j)=reshape(dW(i,:,j),[1 n])*f;
        end
    end
end
%================ Example W in 2 x 2 ======================================
% W = W(x1,x2)
% W = [W11(x1,x2) W12(x1,x2)]
%     [W21(x1,x2) W22(x1,x2)]
% ==================== Jacobian (line3) ===================================
% dW_1 = [dW11,1 dW12,1]
%        [dW21,1 dW22,1]
% dW_2 = [dW11,2 dW12,2]
%        [dW21,2 dW22,2]
% ==================== Shift (line6) [cellarray] ==========================
% dW_1 = [dW11,1 dW11,2] = dW11/dx
%        [dW12,1 dW12,2] = dW12/dx
% dW_2 = [dW21,1 dW21,2] = dW21/dx
%        [dW22,1 dW22,2] = dW22/dx
% ==================== Reshape (line9) [vector] ==========================
% i=1,j=1               i=1,j=2
% Wd(1,1)=dW11 * f      Wd(1,2)=dW12 * f
% i=2,j=1
% W(2,1)=dW21 * f       Wd(2,2)=dW22 * f
