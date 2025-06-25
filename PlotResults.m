%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Exercise for the selection of candidates for the call     %%%
%%%      ‘VAC-2021-42 - PhD Position in CIMNE MARINE’            %%%
%%%                                                             %%%
%%%                   Mohammad Sadegh Eshaghi                   %%%
%%%                                                             %%%
%%%              The function for plotting result               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotResults(t1,y,name)

    figure;
    
    % t and y
    subplot(2,2,1);
    plot(y,'k');
    hold on;
    plot(t1,'r:');
    legend('Outputs','Targets');
    title(name);
    
    % Correlation Plot
    subplot(2,2,2);
    plot(t1,y,'ko');
    hold on;
    xmin=min(min(t1),min(y));
    xmax=max(max(t1),max(y));
    plot([xmin xmax],[xmin xmax],'b','LineWidth',2);
    R=corr(t1',y');
    title(['R = ' num2str(R)]);
    
    % e
    subplot(2,2,3);
    e=t1-y;
    plot(e,'b');
    legend('Error');
    MSE=mean(e.^2);
    RMSE=sqrt(MSE);
    title(['MSE = ' num2str(MSE) ', RMSE = ' num2str(RMSE)]);
    
    subplot(2,2,4);
    histfit(e,50);
    eMean=mean(e);
    eStd=std(e);
    title(['\mu = ' num2str(eMean) ', \sigma = ' num2str(eStd)]);
    
end