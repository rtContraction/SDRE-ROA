function B=SDRE_SDP2Sym(x,A)
% ==================  from SDPVAR to SYMBOLIC  ============================
% IMPORTANT.  IMPUT VECTOR X must contain all n elements
% write the value of A with respect to the variable x
A=str2sym(sdisplay(A));           % variable with  parenthesis inside x(i)
A=string(A); % change to string, preparing to delate ()
h=str2sym(sdisplay(x)); % same as A
h=string(h); %same as A
X = sym('x',[1 length(x)]); % vector with the same variables as x
X=string(X)'; % this element contains only xi = x1,x2...,xn
A=replace(A,h,X); % this string does not contain parenthesis x(i) -> xi
B=str2sym(A); %save the vector as a sym variable
% checked by j.perez 29.10.2020
end