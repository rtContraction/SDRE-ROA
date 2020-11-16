function x1=SDRE_NextState(tau,model,x0,u0)
%% RK45 method
k1=tau*model(0,x0,u0);
k2=tau*model(tau/2,x0+k1/2,u0);
k3=tau*model(tau/2,x0+k2/2,u0);
k4=tau*model(tau,x0+k3,u0);
x1=x0+(k1+2*k2+2*k3+k4)/6;
end