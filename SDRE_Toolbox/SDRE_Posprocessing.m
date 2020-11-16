function [Wsum,v_met_red,Wsumf]=SDRE_Posprocessing(Wnum,v_met,x,tol,ssm)
% =================  Search for non-zero elements  ========================
n = length(x);
Wtest = clean(Wnum,tol); % eliminates elements which coeff < tol
w_met_index = getvariables(Wtest);
v_met_index = getvariables(v_met); % get the *pointers to the variables
h = 0*v_met_index;
for i = 1:length(v_met_index)
    h(i) = getvariables(v_met(i+1)); % first element is always 1
end

[Lia,Locb] = ismember(w_met_index,h); % check Loacb
Locb = Locb+1; %index position + 1 , constant value
if any(any(double(replace(Wtest,x,zeros(n,1)))>tol))
    Locb = [1 Locb];  % Add first element (constant value) and shift array
end

v_met_red = v_met(Locb); % select the elements  included in W(x) 

%sdisplay(v_met)             % original vector with all monomials
%sdisplay(v_met_red)         % reduce version with  val>tol
% =================  Search for non-zero elements  ========================
[EXPW,W_BASE] = getexponentbase(Wtest,v_met_red);
index=full(W_BASE); % get the  coefficients of the matrix after the clean 
l_mon=length(Locb); % get the number of nonlinear tearms 
if size(index,1)>1 % this is true for elements which are matrices
    % =========================  matrix case  =============================
    %%%% Testing %%%
    %sdisplay(Wtest(size(Wtest,1),1)); % element W(N,1)
    %sdisplay(index(size(Wtest,2),:)*v_met_red); % check W(1,N) using index
    %%%% Testing %%%
    Wsum(:,:,l_mon)=zeros(n,n); % nxn for quadratic, change for Y
    for i=1:l_mon
        Wsum(:,:,i)=reshape(index(:,i),n,n); 
        % add code to vectorize as an anonymous (cell?)
    end
    Wsumf=matlabFunction(reshape(sum(index.*SDRE_SDP2Sym(x,v_met(Locb)).'...
    ,2),n,n),'Vars',{SDRE_SDP2Sym(x,x)});  % Function Handle for W
else 
    Wsum=index; %for a simple polynomial element (such as rho)
    %%%% Testing %%%
    %sdisplay(Wtest);
    %sdisplay(index*v_met_red);
    %%%% Testing %%%
    Wsumf=matlabFunction(sum(index.*SDRE_SDP2Sym(x,v_met(Locb)).',2),...
        'Vars',{SDRE_SDP2Sym(x,x)}); % Function Handle for rho
end
v_met_red=matlabFunction(SDRE_SDP2Sym(x,v_met_red),'Vars',...
    {SDRE_SDP2Sym(x,x(ssm))});
end
% checked j.perez 29.10.2020