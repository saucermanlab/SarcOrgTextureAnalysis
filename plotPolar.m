function f = plotPolar(amps)
%% FUNCTION DESCRIPTION
% This function plots the polar graphs as shown in Figure 4c.
%
% INPUTS
% amps : array of sarcomere organization amplitudes corresponding to
% 0 to 180 degrees.
%
% OUTPUTS
% f : function handle

%%
f = figure('units','pixels','position',[0 0 400 400]);

polarplot(deg2rad([0:180 0]),[amps; amps(1); amps(1)],'Color',[1 0 0]); hold on;
p = gca;
p.ThetaLim = [0 180];
p.FontSize = 12;
p.Color = 'none';
p.ThetaTickLabel = {'0°','30°','60°','90°','120°','150°','180°'};
p.ThetaColor = [0 0 0];

p.RLim = [0 0.6];
p.RColor = [0 0 0];
p.GridColor = [0 0 0];
p.GridAlpha = 1;
end