clear;clc
virusType = 'ArchT';
% PCA with spheres
%The input matrix is trials across odors (say 5 trials X 6 odors  = 30) X
%neurons
%% params
numOfOdors = 6;
trialNum = 5;
% rand a color matrix
% clrMat = rand(numOfOdors, 3);
    clrMat = [0.2984    0.9308    0.7855; ...
    0.3330    0.1549    0.7031; ...
    0.7311    0.7012    0.2633; ...
    0.0999    0.8605    0.0948; ...
    0.2058    0.8897    0.9785; ...
    0.7094    0.2452    0.7921];
%% load the data
figure;
for kk = 1:2
    if kk == 2 % load the odor+light matrix
        if strcmpi(virusType, 'ChR2')%ChR2 exps.
            if ismac
                load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/matOdorLight.mat')
            else
                load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\matOdorLight.mat')
            end
        else% ArchT experiments
            if ismac
                load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/matOdorLightArchT.mat')
            else
                load( 'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\matOdorLightArchT.mat')
            end
        end
        
    elseif kk==1 % load the odor matrix
        
        if strcmpi(virusType, 'ChR2')%ChR2 exps.
            if ismac
                load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/matOdor.mat')
            else
                load('E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\matOdor.mat')
            end
        else% ArchT experiments
            if ismac
                load('/Volumes/My Passport/Tal/Experiments/2019/InhibitionOfGC_GAD/dataAnalysis/matOdorArchT.mat')
            else
                load( 'E:\Tal\Experiments\2019\InhibitionOfGC_GAD\dataAnalysis\matOdorArchT.mat')
            end
        end
    end
    
    
    
    %% PCA
    [pc,sc,g]=pca(mat);
    
    subplot(1,2,kk)
    hold all;
    for jj=1:numOfOdors% Run over all odors
        
        % Find the PC values of all current trials
        currValues = [sc((jj-1)*trialNum+1:(jj)*trialNum,1), ...
            sc((jj-1)*trialNum+1:(jj)*trialNum,2),sc((jj-1)*trialNum+1:(jj)*trialNum,3)];
        
        %MEAN+S.E.M
        meanCurrValues = mean(currValues);%In X Y and Z axes.
        seCurrValues = std(currValues)./sqrt(size(currValues,1));%In X Y and Z axes.
       
        
        %Compute the 3d dimensions of the sphere using the SE.
        [X,Y,Z] = sphere;
        X = X*seCurrValues(1);
        Y = Y*seCurrValues(2);
        Z = Z*seCurrValues(3);
        
        %Plot the sphere
        surf(X+meanCurrValues(1),Y+meanCurrValues(2),Z+meanCurrValues(3),  'FaceColor', 'none', 'EdgeColor', clrMat(jj,:))
         zlim([-20 25])
        
        %Scatter all PCA values
        scatter3(currValues(:,1),currValues(:,2), currValues(:,3), 30, clrMat(jj,:), 'fill')
        
        %Clear
        currValues = [];meanCurrValues = [];seCurrValues = [];
        X = []; Y = [];Z = [];
        
    end
    %Plotting parameters
    xlabel('PC1','fontSize',14); ylabel('PC2','fontSize',14); zlabel('PC3','fontSize',14)
    set(gca,'fontSize',14)
    box off;
    set(gca,'tickdir','out','ticklength',get(gca,'ticklength')*2);
    set(gca,'linewidth',.25)
    xlim([-25 40]);ylim([-15 30])
    
    %Titles
    if kk==1
        title('Odor')
    else
         if strcmpi(virusType, 'ChR2')%ChR2 exps.
            title('+ChR2')
         else
             title('ArchT')
         end
    end
    
end
set(gcf, 'Renderer','painters' )