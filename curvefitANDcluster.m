%% Fit Gaussian and use Kmeans clustering to determine groups of cell types encoding interval time

% Load a data set from 'sampledata' folder (one is data using a 300 ms
%   interval, the other is data from presenting an 800 ms interval)
% Data includes 2 graphs and interval-averaged data for each cell, smoothed
% and normalized (i.e., data has been pre-processed)

% Procedures:
% For each neuron, fit a gaussian function to averaged interval data to determine parameters of
%   temporal-associated neural code, with sigma as a free parameter and alpha
%   and beta constrained according to normalization and peak activity window,
%   respectively
% Use fit parameters along with one additional measure (initial neural response to
%   stimulus) as input variables and use unsupervised learning procedure 
%   (Kmeans Clustering) to determine categories of time-coding strategies used by the neurons


%% Initialize some additional variables and user defined options

smoothmore=0;   % consider additional smoothing for noisy data, set to number of passes desired, '0' results in no extra smoothing
                % re-normalizes min-max scaling after smoothing
smoothwindow=5; % number of bins for hanning smoothing window
plotfits=0;     % '1' to plot final fits for each cell, advised if addition smoothin invoked 

kmax=10;        % max number of k clusters to test for cluster algorithm

ISI=[125,200,300,500,800,2000];     % for color-specific graphing based on data set
findcolor=find(ISI==stim_isi,1);
colorcode=[1 0 0;1 0.65 0;0 0 1;0 1 0;1 0 1;0 0 0];
    

%% Fit Gaussian to interval average of each cell
% Constrain alpha to min-max scale (0-1)
% Center fitting window (length = interval length) on time of bin of max response in interval
% Leave sigma unconstrained in positive domain
    
% create output arrays for fit results       
TempTunSigma=zeros(num_cells,1);
TempTunAdjR=zeros(num_cells,1);
TempTunAmplit=zeros(num_cells,1);
TempTunMu=zeros(num_cells,1);
TempTunRMSE=zeros(num_cells,1);
startwin_fits=zeros(num_cells,1);

% create time vector and structure for curvefit function
time=1:binsize:(stim_isi*numcycles);
time=transpose(time);

% create figure if plotting is invoked
if plotfits>0
    plothandle=sprintf('Fig_PlotFits_%d=figure;',stim_isi);
    eval(plothandle);
    set(gcf,'Position',[1    66   749   889]);
    plotcount=1;
end

% loop to pull each cell from data matrix and fit gaussian
fitwin=stim_isi;    % set fit window to length of one interval
for i=1:num_cells
    
    fprintf('Fitting gaussian, neuron #%d\n',i);
    
    % pull activity from one neuron and structure for curvefit function,
    % option to smooth further and re-normalize
    test=iavg_cellmtx(i:i,:);
    if smoothmore>0
        for smoothcount=1:smoothmore
            test=smoothdata(test,'movmean',smoothwindow);
        end
        test=(b-a)*((test-min(test))/(max(test)-min(test)))+a;
    end
    test=transpose(test);
    
    % determine fitting window to center on peak of neuron's activity 
    if binmaxval(i)>=win_shift/2
        startwin=(binmaxval(i)-(win_shift/2))*binsize;       
    else
        startwin=(binmaxval(i)+win_shift-(win_shift/2))*binsize; 
    end

    % use try-catch strategy to deal with errors resulting from failed fits
    try
        % set fit window to exclude data outside that timeframe, centered on bin of maximal response as defined above
            data2excludefit=excludedata(time,test,'domain',[startwin startwin+fitwin-1]);
        % attempt gaussian fit algorithm, constrain alpha to min-max scale and beta (mu) to fit window, leave sigma unconstrained as variable of main interest
            [fitobject,gof]=fit(time,test,'gauss1','Lower', [0 startwin 0],'Upper', [1 startwin+fitwin Inf],'Exclude',data2excludefit);   
        % save final fit results, including start of fit window
            startwin_fits(i,1)=startwin;
            TempTunSigma(i,1)=fitobject.c1;
            TempTunAdjR(i,1)=gof.adjrsquare;
            TempTunAmplit(i,1)=fitobject.a1;
            TempTunRMSE(i,1)=gof.rmse;
            TempTunMu(i,1)=fitobject.b1;
            % if second cycle used to fit mu, subtract so relative to single interval
            if fitobject.b1>stim_isi+responsewindow_4sort
                TempTunMu(i,1)=fitobject.b1-stim_isi;
            else
                TempTunMu(i,1)=fitobject.b1;
            end
    catch
        % if fit fails, record '0' for all output
        TempTunSigma(i,1)=0;
        TempTunAdjR(i,1)=0;
        TempTunAmplit(i,1)=0;
        TempTunMu(i,1)=0;
        TempTunRMSE(i,1)=0;
    end
    
    % plot data with fit if invoked
    if plotfits==1
        subplot(ceil(num_cells/4),4,plotcount);   
            plot(time,test,'-k');
            hold on;
            plot(fitobject,'r');
            s=findobj('type','legend');
            delete(s);
            axis([0 (stim_isi*numcycles) 0 1]);
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            xlabel('')
            ylabel('')
            hold off;
        plotcount=plotcount+1;
    end
    
end

% plot sigma (temporal tuning width) vs mu (time of peak within interval) and perform linear fit
% to quantify scaling of function width over interval, constrain fit to interval
    data2excludefit=excludedata(TempTunMu,TempTunSigma,'range',[1 stim_isi]);
    [fitobject,~]=fit(TempTunMu,TempTunSigma,'poly1','Exclude',data2excludefit);     %TempTunWidth_sorted==0); %TempTunWidth_sorted>stim_isi,
    for y=1:num_cells
        testfit(y,1)=(fitobject.p1*TempTunMu(y,1))+fitobject.p2; %#ok
    end
    m_scale=fitobject.p1;
    fit_SigmaMu=testfit;
    plothandle=sprintf('Fig_muVSsigma_%d=figure;',stim_isi);
    eval(plothandle);
        scatter(TempTunSigma,TempTunMu,'k','filled');
        hold on;
        plot(fit_SigmaMu,TempTunMu,'-r','Color',colorcode(findcolor:findcolor,:),'LineWidth',3);
        set(gca, 'YDir','reverse','TickDir','out')
        set(gcf,'Position',[799   690   440   265])
        axis([-50 stim_isi+75 -50 stim_isi+75]);
        xlabel('Sigma: function width (msec)');
        ylabel('Mu: time to peak of function (msec)');
        title(['ISI ',num2str(stim_isi),': Sigma vs Mu']);
        hold off;
        
 % plot Mu vs measure of fit quality (rmse) to determine if poorer fits are related to time of peak activity during the interval 
    plothandle=sprintf('Fig_pkFuncMuVSrmse_%d=figure;',stim_isi);
    eval(plothandle);
        scatter(TempTunRMSE,TempTunMu,'k','filled');
        set(gca, 'YDir','reverse','TickDir','out')
        set(gcf,'Position',[1241    690         440         265])
        xlim([min(TempTunRMSE)-0.05 max(TempTunRMSE)+0.05]);
        ylim([0 stim_isi+75]);
        xlabel('Fit RMSE: unexplained error');
        ylabel('Mu: time to peak of function (msec)');
        title(['ISI ',num2str(stim_isi),': RMSE vs Mu']);
        hold off

% Calculate last measure of cell response as Difference between first point and min observed within first 100 ms
    min100=min(iavg_cellmtx(:,1:100/binsize),[],2);
    Diff=min100-iavg_cellmtx(:,1:1);
% Plot Scatterplot of Diff vs Mu and bivariate histograpm to visualize 
      plothandle=sprintf('Fig_pkFuncMuVSrespDiff2_%d=figure;',stim_isi);
      eval(plothandle);
          h=scatterhist(Diff,TempTunMu,'Location','SouthEast','Direction','out','Color',colorcode(findcolor:findcolor,:),'LineWidth',2,'Marker','o','NBins',[10,10]);
          axis([-1.05 0.05 0 stim_isi+50]);
          set(h,'TickDir','out','YDir','reverse','XDir','reverse')
          set(gcf,'Position',[777    68   343   297])
          title(['ISI ',num2str(stim_isi),': Response Diff vs Mu']);
          hold off;

      PkbinXDiff=horzcat(TempTunMu,Diff);
      plothandle=sprintf('Fig_DiffVStempMu_%d=figure;',stim_isi);
      eval(plothandle);
      hist3(PkbinXDiff,'FaceColor',colorcode(findcolor:findcolor,:),'FaceAlpha',0.5);
      xlabel('Mu: time to function peak (msec)');
      ylabel('Response Difference (inhib to excit)');
      zlabel('Count');
      title(['ISI ',num2str(stim_isi),': Response Diff vs Mu']);
      set(gcf,'Position',[1121          66         560         420]);
      hold off;
      
%% Use gaussian fit parameters and response Diff measure to cluster cell response types using Kmeans Clustering
% Must rescale each variable, using min-max (0-1) rescaling for this

fprintf('Check for outliers and scale variables.../n');

% Use IQR outer fence method to detect outliers and replace with fence value
% Used only on unconstrained paramter sigma
    [TempTunSigma_outmod]=FindIQRfence(TempTunSigma);
    
% Create data matrix for cluster analysis, each column is variable, each
% row is data for a given neuron
    DataMtx=horzcat(Diff,TempTunSigma_outmod,TempTunMu); 
    numvar=size(DataMtx,2);
    
% Rescale each variable to min-max (0-1)
    for i=1:numvar
        DataMtx_scaled(:,i:i)=(b-a)*((DataMtx(:,i:i)-min(DataMtx(:,i:i)))/(max(DataMtx(:,i:i))-min(DataMtx(:,i:i))))+a;
        checkmin=min(DataMtx_scaled(:,i:i));
        checkmax=max(DataMtx_scaled(:,i:i));
        if checkmin~=a || checkmax~=1
            fprintf('RESCALE VARIABLE FAIL: VAR %d\n',i)
            return
        end
    end
    fprintf('%d variables rescaled: %d - %d\n',i,a,b)
    
% iterate kmeans with incrementing k number of clusters, 100 replicates each
% save output for each run
fprintf('Performing iterative clustering for elbow method of k...\n');
    Dsumsq=NaN(kmax,1);
    D_Idx=NaN(num_cells,1);
    cellID=1:num_cells;
    for i=1:kmax
        [Idx,C,SumD,D]=kmeans(DataMtx_scaled,i,'Replicates',100);
        for j=1:num_cells
            D_Idx(j,1)=D(j,Idx(j));
        end
        Dsumsq(i,1)=sum(D_Idx.^2);
        cmdstr=sprintf("KClust_k%d=struct('ClustID',Idx,'ClustCtr',C,'Distance',D,'DistIdx',D_Idx);",i);
        eval(cmdstr);
    end
    
% plot results to determine best k based on elbow method
    Fig_KsumDsqXk=figure;
        plot(1:kmax,Dsumsq,'-ok');
        set(gcf,'Position',[194   605   499   350])
        xlim([0 kmax]);
        ylim([0 max(Dsumsq)+(max(Dsumsq)*.1)]);
        xlabel('k (number of clusters)');
        ylabel('sum D^2');
        title('Kmeans: k x sum D^2');
        hold off

% query user for k
    [k_elbow]=kQuery;
% Define result dataset based on k
    cmdstr=sprintf('KClust=KClust_k%d',k_elbow);
    eval(cmdstr);

% plot assigned clusters
    ClusterIDplot=figure;
    for i=1:k_elbow
        scatter3(DataMtx(KClust.ClustID==i,1:1),DataMtx(KClust.ClustID==i,2:2),DataMtx(KClust.ClustID==i,3:3),'ok','filled','MarkerFaceColor',color_data(i:i,:));
        hold on;
    end
    xlabel 'Response Difference'
    ylabel 'Temporal Tuning Width (sigma)'
    zlabel 'Temporal Tuning Peak Time (mu)'
    title 'Cluster Assignments and Centroids'
    hold off
 
% Query User regarding plotting different k result
    query=1;
    while query==1
        [query]=UserQueryConfirm;
        if query==1
            [k_elbow]=kQuery;
            cmdstr=sprintf('KClust=KClust_k%d',k_elbow);
            eval(cmdstr);
            % plot new assigned clusters
                ClusterIDplot=figure;
                for i=1:k_elbow
                    scatter3(DataMtx(KClust.ClustID==i,1:1),DataMtx(KClust.ClustID==i,2:2),DataMtx(KClust.ClustID==i,3:3),'ok','filled','MarkerFaceColor',color_data(i:i,:));
                    hold on;
                end
                xlabel 'Response Difference'
                ylabel 'Temporal Tuning Width (sigma)'
                zlabel 'Temporal Tuning Peak Time (mu)'
                title 'Cluster Assignments and Centroids'
                hold off
        elseif query==0
            return;
        end
        
    end
    
    
    
%% Functions

% Find IQR outer fence value and replace outliers, use conservative fence 
function [dataout]=FindIQRfence(data) 

    IQR=quantile(data,0.75)-quantile(data,0.25);
    IQRfence=quantile(data,0.75)+(IQR*4);
    
    data(data>IQRfence,1)=IQRfence;
    dataout=data;
    
end

% Query User for k
function [k_elbow]=kQuery

    prompt={'Enter k (ideal number of clusters based on plot elbow)'};
    dlgtitle='Choose k';
    definput={'3'};
    answer=inputdlg(prompt,dlgtitle,[1 40],definput);
    k_elbow=str2double(answer);
    
end


% Confirm k query
function [query]=UserQueryConfirm
    
    answer = questdlg('Would you like to plot a different k result?','Replot Request','No, looks good!','Yes, please!','No, looks good!');
    % Handle response
    switch answer
        case 'No, looks good!'
            query=0;

        case 'Yes, please!'
            query = 1;
    end
    
end