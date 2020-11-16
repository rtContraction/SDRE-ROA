function output=SDRE_SaveTikz(system,intLength,Platform)
output=1;
clc
close all

if strcmp(system,'secondOrder')
    % Plots for the second order system
    % The labels are automatically read in Latex
    xLength=5;yLength=5;% 
    finalTime=3; % less or equal to sdre.simT
    lineWidth=1.2;
    figNames={'states.fig','rate.fig','norm.fig','region.fig'};
    for i=1:length(figNames)
    openfig(figNames{i});
    cleanfigure; % important to export to tikz
    end
    %%
    empty={' '}; % label for all grey plots
    strCopy=repmat(empty,1,intLength); 
    strColor={'$\mbx_1$','$\mbx_2$','$\mbx_3$','$\mbx_4$','$\mbx_5$'};
    figure(2);
    figLabel=flip(findall(gcf,'Type','Line')); %the order is reversed
    strFix={'$\cRate_C$','$\cRate$'};
    legendstr=[strFix,strCopy,strColor];
    selectLabel=~strcmp(empty,legendstr);
    legend(figLabel(selectLabel),...
        legendstr{selectLabel},'interpreter','latex')
    legend('Location','northeast')
    ylabel('$\cRate(\mbx)$','interpreter','latex')
    xlabel('time','interpreter','latex')
    axis([0 finalTime 0 35])
    grid('on')
    box on
    SDRE_ExportTikz('ContractionRate.tikz',Platform)

    figure(3);
    figLabel=flip(findall(gcf,'Type','Line'));
    legendstr=[strCopy,strColor];
    selectLabel=~strcmp(empty,legendstr);
    legend(figLabel(selectLabel),...
        legendstr{selectLabel},'interpreter','latex')
    legend('Location','northeast')
    xlabel('time','interpreter','latex')
    ylabel('$\parallel \mbx \parallel_{\mbM}$','interpreter','latex')
    grid('on')
    axis([0 finalTime 0 40])
    box on
    SDRE_ExportTikz('Convergence.tikz',Platform)

    figure(4);
    figLabel=[flip(findall(gcf,'Type','Line'));...
        flip(findall(gcf,'Type','Contour'))];
    strFix={'$\ROABracci$','$\ROAChang$'};
    strCopy=repmat(empty,1,intLength+2);
    %strColor={'$\mbx_1$',' ','$\mbx_2$',' ','$\mbx_3$',' ',...
    %'$\mbx_4$',' ','$\mbx_5$',' ','$\ROAZero$','$\ROAOpt$'}; % trajectories
    strColor={' ',' ',' ',' ',' ',...
    ' ',' ',' ','$\ROAZero$','$\ROAOpt$'}; %only ROA %forget plot  !
    legendstr=[strFix,strCopy,strColor];
    selectLabel=~strcmp(empty,legendstr);
    axis([-xLength xLength -yLength yLength])
    yticks(-xLength+1:2:xLength-1)
    xticks(-yLength+1:2:yLength-1)
    box on
    legend(figLabel(selectLabel),...
        legendstr{selectLabel},'interpreter','latex')
    legend('Location','southwest')
    xlabel('$x_1$','interpreter','latex')
    ylabel('$x_2$','interpreter','latex')
    %SDRE_ExportTikz('RegionAttraction.tikz',Platform)
    
    %%
    close(figure(4))
    figure(4);
    load('sdreTemp.mat','sdre','sys')
    H=[];
    xPoints=20;yPoints=20; % set S where the pointwise properties are checked
    [X,~]=meshgrid(linspace(-xLength,xLength,xPoints+1),...
        linspace(-yLength,yLength,yPoints+1));
    for i=1:size(X,1)
       for j=1:size(X,2)
           H=[H;[X(1,i),X(1,j)]];
       end   
    end
    cShade=[102,255,255]/255;
    plot(0,0,':m','linewidth',lineWidth)
    hold on
    plot(H(:,1),H(:,2),'Marker','x','Color',cShade,'LineStyle','none')
    bracciData([xLength;yLength],lineWidth);
    sdre.x0=[2;2];
    color=6;
    SDRE_Simulation(sdre,sys,color,0)
    figLabel=[flip(findall(gcf,'Type','Line'));...
        flip(findall(gcf,'Type','Contour'))];
    strFix={'$\ROABracci$','$\mbx_g \in \ROAChang$'};
    strCopy=repmat(empty,1,length(figLabel)-2);
    legendstr=[strFix,strCopy];
    selectLabel=~strcmp(empty,legendstr);
    legend(figLabel(selectLabel),...
        legendstr{selectLabel},'interpreter','latex')
    legend('Location','southwest')
    xlabel('$x_1$','interpreter','latex')
    ylabel('$x_2$','interpreter','latex')
    axis([-xLength xLength -yLength yLength])
    yticks(-xLength+1:2:xLength-1)
    xticks(-yLength+1:2:yLength-1)
    box on
    cleanfigure;
    %SDRE_ExportTikz('BracciChang.tikz',Platform)
elseif strcmp(system,'TWIP')
    figNames={'TWIPMetric.fig','TWIPNormPosition.fig'...
                            ,'TWIPNormVelocity.fig','TWIPNormInput.fig'};
    for i=1:length(figNames)
    openfig(figNames{i});
    cleanfigure;
    end
    finalTime=1;
    
    figure(1); 
    ylabel('rate,norm','interpreter','latex')
    legend({'$\cRate$','$\cRate(\mbx)$',...
        '$ \parallel \mbx \parallel_{\mbM}$'},'interpreter','latex')
    xlabel('time','interpreter','latex')
    %axis([0 finalTime 0 30])
    grid('on')
    SDRE_ExportTikz('WIPMetric.tikz',Platform)
    
    figure(2); 
    legendstr={'\nameLQRA','\nameLQRB','SDRE'};
    legend(legendstr,'interpreter','latex')
    xlabel('time','interpreter','latex')
    ylabel('$\parallel \posV \parallel$','interpreter','latex')
    grid('on')
    %axis([0 finalTime 0 30])
    SDRE_ExportTikz('WIPNormPosition.tikz',Platform)
    
 	figure(3); 
    legend(legendstr,'interpreter','latex')
    xlabel('time','interpreter','latex')
    ylabel('$\parallel \velV \parallel$','interpreter','latex')
    grid('on')
    %axis([0 finalTime 0 30])
    SDRE_ExportTikz('WIPNormVelocity.tikz',Platform)
    
    figure(4);
    legend(legendstr,'interpreter','latex')
    xlabel('time','interpreter','latex')
    ylabel('$\parallel \mbu \parallel$','interpreter','latex')
    grid('on')
    axis([0 finalTime 0 Inf])
    SDRE_ExportTikz('WIPNormControl.tikz',Platform)
end
end