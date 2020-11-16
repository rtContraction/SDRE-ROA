%% ROA Estimate Bracci 

xPoints=50;yPoints=50; % set S where the pointwise properties are checked
xLength=5;yLength=5;% 8 Metric for control
[PclA,lyapFun]=bracciData([xLength;yLength],1.2,[xPoints;yPoints]);
