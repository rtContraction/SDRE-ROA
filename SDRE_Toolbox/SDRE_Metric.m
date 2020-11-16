function [W,C,v]=SDRE_Metric(x,sys)
% ==============================  Initialization ==========================
n=length(x);
k=sys.met_order;
index=sys.v_ssm;
const=sys.v_ssm_const;
o=length(const); % number of states that will be constant
l=o*(o+1)/2; % number of constant entries
r=n*(n+1)/2; % number of entries requeried for each  symmetric Matrix
h=nchoosek(length(index)+k,k); % number of combinations to generate h monom
                               % (variables in x + degree), degree
P=sdpvar(r-l,1); % P_i is a scalar sdpvar
Q=sdpvar(l,1);   % sdpvars for constant entries of W
C=sdpvar(h,r-l); % C_i is a vector of sdpvar
for i=1:(r-l)      
    [P(i,1),C(:,i),v] = polynomial(x(index),k); % pol with monom of x(m)
end
W=sdpvar(n,n,'symmetric');
j=1;q=1;p=1; %column, row index for generating a Symmetric Matrix 
% while(j<=n) % W11 -> W12 M22 --> W13 W23 W33.. (First columns, add (i=1))
%     %i
% while( i<=j)
%     %j
%     M(i,j)=P(q,1);
%     i=i+1;
%     q=q+1;
% end
% j=j+1;
% i=1;
% end   
% ===========================  Parametric W ===============================
while(j<=n)     % W11 W12 W13-->  W22 W23  -> W33... First rows
    i=j;
while(i<=n)
    %j
    if ismember(i,const) && ismember(j,const)
        W(j,i)=Q(q,1);
        q=q+1;
    else
        W(j,i)=P(p,1);
        p=p+1;
    end
    i=i+1;
end
j=j+1;
end
W= triu(W,0) + triu(W,1).'; % W=W' using the Upper triangular From
if isempty(C)
    v=1;
else
    C = reshape(C,[length(C(:,1))*length(C(1,:)),1]); %  from matrix to vector
end
C = [C;Q];  %add variables for const entries of matrix to the vector
% checked by j.perez 29.10.2020

