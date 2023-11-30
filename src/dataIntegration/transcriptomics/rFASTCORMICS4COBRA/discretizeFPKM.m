function [discretized] = discretizeFPKM(fpkm,colnames,figflag, p)
%The function discretizes the data by fitting gaussian distribution 
%please check Pacheco et al for more details

% USAGE:
%
%   [discretized] = discretizeFPKM(fpkm,colnames,figflag, p))

% INPUTS:
% fpkm:                  fpkm values for the samples, size(fpkm,1) =
%                        lenghth(geneIDs), size(fpkm,2)= length(colnames)
% colnames:              cell array with the sample names


% OPTIONAL INPUTS:

% figflag                print and save Figures if 1 do not print if 0
% path                    Folder in which the plots should be saved



% OUTPUTS:

% discretized:           discretized values in 1 (expression), -1 (no
%                        expression, 0 unknown expression status

% .. Authors:
%       - Maria Pires Pacheco, Thomas Sauter, 2016, University of Luxembourg
%       - Maria Pires Pacheco, Thomas Sauter, 2022, adaptation of the code to the Cobra toolbox


if nargin <4
      p=genpath(pwd);
      p=strcat(p,'\');
      
end
      q=strrep(p,'\','/');

if nargin <3
    figflag = 0;
end

if istable(fpkm)
    fpkm = table2array(fpkm);
end

signal = log2(fpkm); % log2-transform the data

signal(isinf(signal)) = -1e6;
discretizedKeep= zeros(size(fpkm,2),size(fpkm,1));

for j=1:size(fpkm,2) % Discretize the data by creating half-gaussians

    signalSample = signal(:,j);
    data_keep = fpkm(:,j);
    signalSample = signalSample(signalSample > -1e6);
    
    if figflag; figure  %create folder to save Figures if doesnt exist        
        if ~exist([q,'Figures/Discretization/'], 'dir')
            mkdir([q,'Figures/Discretization/'])
        end
        
        
        hist((signalSample),100);
        ylabel('abundance');xlabel('log2(FPKM)');title({'Histogram of abundance',['Sample: ', colnames{j}]},'Interpreter','none');
        saveas(gcf,[p,'Figures\Discretization\histogram_', colnames{j},'.png'])
    end
    
    [probabilityEstimate,xi] = ksdensity((signalSample)); % Density plot

    if figflag;figure; hold on
        plot(xi,probabilityEstimate,'-k','LineWidth',2);
        title({'Density curve for sample: ', colnames{j}},'Interpreter','none');
        ylabel('Density'); xlabel('log2(FPKM)');
        title({'Density plot',['Sample: ', colnames{j}]},'Interpreter','none');
        legend({'Signal'},'Location','best')
        saveas(gcf,[p,'Figures\Discretization\density_', colnames{j},'.png']);
    end
    
    peak_max_idx = (find(probabilityEstimate == max(probabilityEstimate))); % find the maximum of the density curve
    maxRside = probabilityEstimate(peak_max_idx:end); % copy right-most side of density curve
    maxLside = flip(maxRside); % mirror the right side
    hybrid = [maxLside(1:end-1), maxRside];
    
    
    hybrid_curve = zeros(numel(probabilityEstimate),1);
    if numel(hybrid)> numel(probabilityEstimate)
        hybrid_curve(end+1-numel(hybrid(end-100+1:end)):end) = hybrid(end-100+1:end);
        x = xi;
    else
        hybrid_curve(end-numel(hybrid)+1:end) = hybrid;
        x = xi;
    end
    
    if figflag
        plot(x,hybrid_curve,'b--','LineWidth',2);
        legend({'Signal','hybrid curve right'},'Location','best')
        title({'Density plot with hybrid curves',['Sample: ', colnames{j}]},'Interpreter','none');
    end
    
    clear peak_max_idx maxLside maxRside hybrid
    
    rest_curve = probabilityEstimate - hybrid_curve';   % create left curve (main curve - rightmost curve)
    rest_curve(rest_curve<0.0001)=0;
    
    if figflag; plot(xi,rest_curve, 'r--','LineWidth',2);
        legend({'Signal','hybrid curve right','hybrid curve left'},'Location','best')
        saveas(gcf,[p,'Figures\Discretization\density_maxfit_', colnames{j},'.png']);
    end
    
    [fitresultR, ~] = createFitR(xi, hybrid_curve); % fit the curves
    [fitresultL, ~] = createFitL(xi, rest_curve);
    
    if figflag
        plot(fitresultR,'b');
        plot(fitresultL),'g';
        title({'Density plot with fitted curves',['Sample: ', colnames{j}]},'Interpreter','none');
        ylabel('Density'); xlabel('log2(FPKM)');
        legend({'Signal', 'hybrid curve right', 'hybrid curve left','fitted curve right','fitted curve left'},'Location','best')
        saveas(gcf,[p,'Figures\Discretization\density_fitted_', colnames{j},'.png']);
        close all;
    end
    
    sigma1 = fitresultR.c1/sqrt(2);    %% zFPKM transform the data and plot the density
    mu1 = fitresultR.b1; % x-value of right max
    
    zFPKM = (signalSample-mu1)/sigma1;
    [yFPKM,xFPKM] = ksdensity(zFPKM);
    

   
    zFPKM = (signal(:,j)-mu1)/sigma1; %transform sample into zFPKM
    
    %sigma2  = fitresultL.c1/sqrt(2);
    mu2     = fitresultL.b1;
    mu2_t   = (mu2-mu1)/sigma1; % x-value of left max
    if figflag
        plot([mu1, mu1],[0, max(probabilityEstimate)])
        plot([mu2, mu2],[0, max(probabilityEstimate)])
        legend({'Signal', 'hybrid curve right', 'hybrid curve left',...
            'fitted curve right','fitted curve left',...
            ['expression threshold = ',num2str(mu1)],['inexpression threshold = ',num2str(mu2)]},'Location','best')
    end
    
    discretized = zFPKM;
    
    e = 0;
    ue = max(-3, mu2_t);
    
    
    if figflag; figure;hold on;
        plot(xFPKM,yFPKM,'-k','LineWidth',2);
        xlabel('zFPKM');
        ylabel('density');
        line([0 0], [0.0 0.3],'color', [0,1,0])
        line([ue ue], [0.0 0.3],'color', [0,1,1]);
        legend({'zFPKM data','expression threshold','inexpression threshold'},'Location','best')
        title({'Expression thresholds',['Sample: ', colnames{j}],['inexpression threshold = ', num2str(ue)]},'Interpreter','none');
        
        saveas(gcf,[p,'Figures\Discretization\zFPKM_', colnames{j},'.png']);
        close all;
    end
    
    
    
    expThreshold   = e;
    unexpThreshold = ue;
    
    expressed     = discretized >= expThreshold;
    not_exp = discretized <= unexpThreshold;
    unknown = discretized < expThreshold & discretized>unexpThreshold;
    
    discretized(unknown)    = 0;
    discretized(not_exp)    = -1;
    discretized(expressed)  = 1;
    
    discretized=(reshape(discretized,size(data_keep,1),size(data_keep,2)));
    clear fitresultL fitresultR data_keep e exp expThreshold mu1 mu2 mu2_t not_exp sigma1 sigma2 signal2 signalSample ue unexpThreshold unknown zFPKM
    discretizedKeep(j,:)=discretized';
end
discretized = discretizedKeep';
clear discretizedKeep j

if figflag;figure; hold on %% Density plot for all samples
    for j=1:length(colnames) 
        signalSample = signal(:,j);
        signalSample = signalSample(exp(signalSample) > 1e-6);
        
        [probabilityEstimate,xi] = ksdensity((signalSample));         % Density plot

        plot(xi,probabilityEstimate,'LineWidth',1);
    end
    title({'Density curves for all samples'});
    ylabel('Density'); xlabel('log2(FPKM)');
    saveas(gcf,[p,'Figures\Discretization\density_all.png']);close all
end

if figflag;figure; hold on % Subplots for all samples
    for j=1:length(colnames) 
        signalSample = signal(:,j);
        signalSample = signalSample(exp(signalSample) > 1e-6);        
        [probabilityEstimate,xi] = ksdensity((signalSample)); % Density plot

        plot(xi,probabilityEstimate,'LineWidth',1);
        title(colnames{j}, 'interpreter','none');
    end
    saveas(gcf, [p,'Figures\Discretization\density_subplots.png']);close all
end
clear j probabilityEstimate signalSample xi
end
