function [Treal, Tboot, allMeans, SE, CI, pValue] = ...
    hBS_Munc13_function(fileName, hFlag, nBoot, myAlpha, pFlag)

% hBS_Munc13_function.m: hierarchical bootstrap of Kaeser lab data
%
% [Treal, Tboot, allMeans, SE, CI, pVal] = 
%       hBS_Munc13_function(fileName, hFlag, nBoot, myAlpha, pFlag)
%
% ex. [Treal, Tboot, allMeans, SE, CI, pVal] = hBS_Munc13_function('dataset 1-sucrose IPSC-v4.xlsx', 0, 10000,0.05, 1);
%
% Inputs:
% - fileName, name of an excel data file
% - hFlag, 1 = hierarchical bootstrap (default); 0 = non-hierarchical bootstrap
% - nBoot, # of bootstrap replicates (default = 1000)
% - myAlpha, significance level for confidence intervals (default = 0.05)
% - pFlag, 1 = display histogram (default = 0, no histogram)
%
% Outputs:
% - Treal, experimental value of test statistic, T
% - Tboot, bootstrap replicates of the test statsitic, T
% - allMeans, grand mean of each group for each bootstrap iteration (4 x nBoot)
% - CI, confidence interval corresponding to myAlpha
% - pValue, probability of H0 given our data
%
% 29 September 2022, RTB wrote it; plane ride to MKE, mom's hip surgery
% 08 October 2022, converted to function & added option for
%    non-hierarchical bootstrap

% When the Kaeser lab knocks out two genes involved in synaptic
% release (RIMs & ELKS), synaptic transmission is reduced by about 80%. Now
% they want to ask whether KO'ing a 3rd gene, Munc13, further reduces
% synaptic transmission. However, since the double KO (RIMs/ELKS) is in a
% different mouse strain from the triple KO (RIMs/ELKS/Munc13), the
% synaptic transmission in each strain's KO must be compared to its own
% control.
% 
% S_1 = Strain 1: RIMs/ELKS-Cre
% S_2 = Strain 2: RIMs/ELKS/Munc13-Cre
% 
% To perform an experiment in a given strain, several mice are killed, the
% hippocampi (HC) are dissected, cells dissociated and pooled together in
% one big primary culture. This culture is then divided into two groups:
% one is treated with a lentivirus containing Cre (= KO); the other is
% treated with a lentivirus containing a dead Cre, 'delta-Cre' (= Control).
% From each culture, multiple cells are tested. Each cell is patch clamped
% and synaptic transmission is tested by measuring the size of the
% post-synaptic current (EPSC) evoked by an action potential. This
% measurment is generally repeated 5 times for each cell.
% 
% C_1 = Condition 1: KO
% C_2 = Condition 2: Control
% 
% Thus each group of 5 related measurements is uniquely identified by FOUR
% numbers, which will be variables (columns) in the Excel file. The Excel
% file will be read in as a Table variable, 'ds':
% 
% ds.Strain [1:2]
% ds.Cond [1:2]
% ds.Batch [1:nBatches], may vary for different experiments
% ds.Cell [1:nCells], may vary for different batches
% 
% For ease of discussion, we'll think of four groups:
% 
% Group A: S_1, C_1 (RIMS/ELKS KO)
% Group B: S_1, C_2 (RIMS/ELKS Control)
% Group C: S_2, C_1 (RIMS/ELKS/Munc13 KO)
% Group D: S_2, C_2 (RIMS/ELKS/Munc13 Control)
% 
% The scientific question is whether knocking out Munc13 (triple mutant,
% S2) causes a greater relative decrease in synaptic transmission as
% compared to the double mutant (S1). So we will define our test statistic,
% T, as:
% 
% T = (mean(Group_A) / mean(Group_B)) / . . . 
%     (mean(Group_C) / mean(Group_D)) 
% 
% The null value for this statistic is 1; our alternate hypothesis is:
%       T > 1
% 
% In this hierarchical bootstrap, we directly estimate the sampling
% distribution of our test statistic by re-sampling, with replacement, from
% the raw data, while preserving the hierarchical relationships in the
% data. For each bootstrap iteration, Group identity (fixed effect) is
% preserved. We do our re-sampling at three nested levels of what would be
% random effects in a linear mixed effects model: "batch," "cell" and
% "sweep" (technical replicates).

%% Read the Excel file into a data table

if nargin < 5, pFlag = 0; end
if nargin < 4, myAlpha = 0.05; end
if nargin < 3, nBoot = 1000; end
if nargin < 2, hFlag = 1; end

if nargin < 1
    error('Error: A data file name must be provided.');
end

ds = readtable(fileName);

% Check the column names
varNames = ds.Properties.VariableNames;

%% Calculate the actual value of our test statistic, T

dsGrpA = ds((ds.Strain == 1) & (ds.Condition == 1),:);  % double KO Cre
dsGrpB = ds((ds.Strain == 1) & (ds.Condition == 2),:);  % double KO control
dsGrpC = ds((ds.Strain == 2) & (ds.Condition == 1),:);  % triple KO Cre
dsGrpD = ds((ds.Strain == 2) & (ds.Condition == 2),:);  % triple KO control

Treal = (mean(dsGrpA.PSC,'omitnan') / mean(dsGrpB.PSC,'omitnan')) / ...
        (mean(dsGrpC.PSC,'omitnan') / mean(dsGrpD.PSC,'omitnan'));

%% Hierarchical bootstrap

% For each group, we first resample (with replacement) from the batches.
% Then, for each of these batches, we resample from the cells; then for
% each cell we resample from the (usually 5) technical replicates. For the
% "sucrose" experiments, there is only one data point per cell.

% Many 'for' loops, and a lot of subsetting of the data (which can probably
% be avoided with more thought), so this might be kinda slow. Let's time it
% for different values of nBoot
tic;
% For nBoot = 100,000, run time of about 45 min.
% For nBoot = 10,000 the run time was around 353 seconds (almost 6 minutes)
% For nBoot = 1,000, run time is 28 seconds
% For nBoot = 100, run time is about 3 seconds

% define constants:
nStrains = 2;   % S_1 = RIMs/ELKS, S_2 = RIMs/ELKS/Munc13
nConds = 2;     % C_1 = KO, C_2 = control

% variable to store the grand mean value of the EPSCs for each group (rows)
% and each bootstrap iteration (columns):
allMeans = zeros(nStrains*nConds, nBoot);

for thisBoot = 1:nBoot
    % Two strains: RIMS/ELKS-Cre and RIMS/ELKS/Munc13-Cre
    for thisStrain = 1: nStrains
        % Two conditions: Cre and control
        for thisCond = 1:nConds
            
            % Temporary variable to hold the resampled EPSC values:
            dataStar = [];
            
            % Grab the subset of the data corresponding to this group:
            dsGrp = ds((ds.Strain == thisStrain) & (ds.Condition == thisCond),:);
            
            % hFlag == 1, do the Hierarchical Bootstrap
            if hFlag
                
                % How many batches of neurons for this group?
                nBatches = length(unique(dsGrp.Batch,'rows'));
                
                % Re-sample with replacement from batches (Note that this
                % assumes batches go from 1:nBatches with no gaps.)
                bStarNdx = unidrnd(nBatches, nBatches, 1);
                
                for thisBatch = 1:nBatches
                    % Grab the subset of the data corresponding to this batch:
                    dsBatch = dsGrp(dsGrp.Batch == bStarNdx(thisBatch),:);
                    
                    % How many cells in this batch?
                    nCells = length(unique(dsBatch.Cell,'rows'));
                    
                    % Re-sample with replacement from cells
                    cStarNdx = unidrnd(nCells, nCells, 1);
                    
                    for thisCell = 1:nCells
                        % dsCell contains all the measurements for one cell:
                        dsCell = dsBatch(dsBatch.Cell == cStarNdx(thisCell),:);
                        
                        % How many technical replicates for this cell?
                        nSweeps = height(dsCell);
                        
                        % For sucrose experiments, there is only one value per
                        % cell
                        if nSweeps == 1
                            dataStar = [dataStar; dsCell.PSC];
                        else
                            % Re-sample with replacement from replicates
                            swStarNdx = unidrnd(nSweeps, nSweeps, 1);
                            dataStar = [dataStar; dsCell.PSC(swStarNdx)];
                        end
                    end
                end
            
            % Non-hierarchical bootstrap
            else
                % resample with replacement from all recordings
                nGrp = height(dsGrp);
                dataStar = dsGrp.PSC(unidrnd(nGrp,nGrp,1));
            end
            
            % Store mean of dataStar in allMeans
            allMeans(((thisStrain - 1)*2) + thisCond, thisBoot) = mean(dataStar,'omitnan');
            
        end
    end
end

% Calculate our test statistic, T, for each bootstrap iteration:
Tboot = (allMeans(1,:) ./ allMeans(2,:)) ./ (allMeans(3,:) ./ allMeans(4,:));

elapsedTimeInSeconds = toc;
disp(['Run time was ' num2str((elapsedTimeInSeconds/60),2) ' min.']);

%% Calculate standard error, confidence intervals and a p-value

% calculate the standard error:
SE = std(Tboot);

% calculate 95% CI using percentile method:
sortedTboot = sort(Tboot);
idxHi = ceil(nBoot * (1 - myAlpha/2));
idxLo = floor(nBoot * (myAlpha/2));
CIhi = sortedTboot(idxHi);
CIlo = sortedTboot(idxLo);
CI = [CIlo, CIhi];

% calculate p-value:
pValue = sum(Tboot <= 1) / nBoot;

if pValue == 0
    pValue = 1 / (nBoot + 1);
end

%% Plot a histogram of our bootstrap distribution

if pFlag
    
    if hFlag
        figure('Name','Hierarchical Bootstrap');
    else
        figure('Name','Non-hierarchical Booststrap');
    end
    %histogram(Tboot);
    hist(Tboot,50);
    hold on
    
    % add a line for our experimental value
    ax = axis;
    h1 = line([Treal, Treal],[ax(3), ax(4)]);
    set(h1, 'Color', [0.6350, 0.0780, 0.1840],'LineWidth',2);
    
    % add dashed lines for confidence intervals
    h2 = line([CIhi, CIhi],[ax(3), ax(4)]);
    set(h2, 'Color', [0.6350, 0.0780, 0.1840],'LineStyle','--','LineWidth',1.5);
    h3 = line([CIlo, CIlo],[ax(3), ax(4)]);
    set(h3, 'Color', [0.6350, 0.0780, 0.1840],'LineStyle','--','LineWidth',1.5);
    
    % label axes, etc.
    xlabel('T^{*}');
    ylabel('Frequency');
    title([fileName '     p(H0|Data) = ' num2str(pValue,3)]);
    ciStr = [num2str((1-myAlpha)*100) '% CI'];
    legend([h1,h2],{'Experimental value',ciStr});
    set(gca,'FontSize',12,'TickDir','out');
    
end