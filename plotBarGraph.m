function plotBarGraph(bad,good)
%% FUNCTION DESCRIPTION
% This function plots bars used in Figure 5. This only plots the actual
% bars and standard deviations.
%
% INPUTS
% bad : array for cells with disorganized sarcomeres
% good : array for cells with organized sarcomeres

figure('units','pixels','position',[50 50 350 800]); hold on;

bar(0.9,mean(bad),1,'FaceColor',[189,189,189]./255,'EdgeColor','none','LineWidth',2);
bar(2.1,mean(good),1,'FaceColor',[99,99,99]./255,'EdgeColor','none','LineWidth',2);

plotstd(0.9,mean(bad),std(bad),0.15,2,0);
plotstd(2.1,mean(good),std(good),0.15,2,0);

xlim([0 3]);
ylabel('Score');

ylim([0 1.75*mean(good)]);

set(gca,'FontSize',24,'Color','none','TickDir','out','LineWidth',2);
set(gca,'XTick',[]);

set(gca,'Position',[0.3 0.05 0.6 0.85]);
end