clear;
min_frequency=110;
min_f=num2str(min_frequency);
Parameters='Three';
day='330';
cd(['C:\Users\Tom\OneDrive\ASU_Thesis\Beams\EDGES\high_band\day' day]); 

%%
% Case 1
% Beta with CMB fit and with CMB adjusted beam correction
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P0_CMB_fit_with_beam_correction_with_cmb_adjust.mat'])
BetaCMBfitP0=Expression1;
LST=BetaCMBfitP0(:,1);
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P5_CMB_fit_with_beam_correction_with_cmb_adjust.mat'])
BetaCMBfitP5=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P180_CMB_fit_with_beam_correction_with_cmb_adjust.mat'])
BetaCMBfitP180=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P185_CMB_fit_with_beam_correction_with_cmb_adjust.mat'])
BetaCMBfitP185=Expression1;

figure(1);clf(1);plot(LST,BetaCMBfitP0(:,2));
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.01 .35 .30 .35])
hold on
plot(LST,BetaCMBfitP5(:,2))
plot(LST,BetaCMBfitP180(:,2))
plot(LST,BetaCMBfitP185(:,2)); ylim([2.35 2.65]); xlim([0 24]); ax=gca; ax.XTick=[(0:4:24)];
grid; ylabel('Beta');xlabel('LST (h)'); 
legend('P0','P5','P180','P185','Location','southeast');
title('With CMB Removed in the Fit and CMB Considered in Beam Adjustment');
text(1,2.4,{['Min F = ' min_f ' MHz'], [Parameters ' Params'], ['Day = ' day]})

set(gcf, 'PaperSize', [7.5 6]); %Set the paper to have width 7.5 and height 6.
set(gcf, 'PaperPosition', [0 0 7.5 6]); %Position plot at left hand corner with width 7.5 and height 6.
savefiledir='C:\Users\Tom\OneDrive\ASU_Thesis\Edges\LOCO Memos\Spectral Index Reports 2015-11 and 2016-06\Figures\';
saveas(gcf, [savefiledir 'FullCMB_PVariations_Fmin' min_f 'MHz_Day' day '_' Parameters '_Params.jpg'], 'jpg');

%%
% Case 2
% Beta without CMB fit and with CMB adjusted beam correction
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P0_no_CMB_fit_with_beam_correction_with_cmb_adjust.mat'])
BetaNoCMBfitP0=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P5_no_CMB_fit_with_beam_correction_with_cmb_adjust.mat'])
BetaNoCMBfitP5=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P180_no_CMB_fit_with_beam_correction_with_cmb_adjust.mat'])
BetaNoCMBfitP180=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P185_no_CMB_fit_with_beam_correction_with_cmb_adjust.mat'])
BetaNoCMBfitP185=Expression1;

figure(2);clf(2);plot(LST,BetaNoCMBfitP0(:,2))
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.28 .35 .30 .35])
hold on
plot(LST,BetaNoCMBfitP5(:,2))
plot(LST,BetaNoCMBfitP180(:,2))
plot(LST,BetaNoCMBfitP185(:,2)); ylim([2.35 2.65]);xlim([0 24]); ax=gca; ax.XTick=[(0:4:24)];
grid; ylabel('Beta');xlabel('LST (h)');
legend('P0','P5','P180','P185','Location','southeast')
title('Without CMB Removed in the Fit and CMB Considered in Beam Adjustment');
text(1,2.4,{['Min F = ' min_f ' MHz'], [Parameters ' Params'], ['Day = ' day]})

%%
% Case 3
% Beta with CMB fit and without CMB adjusted beam correction
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P0_CMB_fit_with_beam_correction_without_cmb_adjust.mat'])
BetaCMBfitNoCMBBeamAdjustP0=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P5_CMB_fit_with_beam_correction_without_cmb_adjust.mat'])
BetaCMBfitNoCMBBeamAdjustP5=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P180_CMB_fit_with_beam_correction_without_cmb_adjust.mat'])
BetaCMBfitNoCMBBeamAdjustP180=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P185_CMB_fit_with_beam_correction_without_cmb_adjust.mat'])
BetaCMBfitNoCMBBeamAdjustP185=Expression1;

figure(3);clf(3);plot(LST,BetaCMBfitNoCMBBeamAdjustP0(:,2));
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.01 .02 .30 .35]);
hold on
plot(LST,BetaCMBfitNoCMBBeamAdjustP5(:,2))
plot(LST,BetaCMBfitNoCMBBeamAdjustP180(:,2))
plot(LST,BetaCMBfitNoCMBBeamAdjustP185(:,2)); ylim([2.35 2.65]);xlim([0 24]); ax=gca; ax.XTick=[(0:4:24)];
grid; ylabel('Beta');xlabel('LST (h)');
legend('P0','P5','P180','P185','Location','southeast');
title('With CMB Removed in the Fit and CMB Not Considered in Beam Adjustment');
text(1,2.4,{['Min F = ' min_f ' MHz'], [Parameters ' Params'], ['Day = ' day]})

%%
% Case 4
% Beta without CMB fit and without CMB adjusted beam correction
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P0_no_CMB_fit_with_beam_correction_without_cmb_adjust.mat'])
BetaNoCMBfitNoCMBBeamAdjustP0=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P5_no_CMB_fit_with_beam_correction_without_cmb_adjust.mat'])
BetaNoCMBfitNoCMBBeamAdjustP5=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P180_no_CMB_fit_with_beam_correction_without_cmb_adjust.mat'])
BetaNoCMBfitNoCMBBeamAdjustP180=Expression1;
load(['Beta_day' day '_Blade_' Parameters '_Params_' min_f '-190MHz_day' day '_P185_no_CMB_fit_with_beam_correction_without_cmb_adjust.mat'])
BetaNoCMBfitNoCMBBeamAdjustP185=Expression1;

figure(4);clf(4);plot(LST,BetaNoCMBfitNoCMBBeamAdjustP0(:,2))
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.28 .02 .30 .35])
hold on
plot(LST,BetaNoCMBfitNoCMBBeamAdjustP5(:,2))
plot(LST,BetaNoCMBfitNoCMBBeamAdjustP180(:,2))
plot(LST,BetaNoCMBfitNoCMBBeamAdjustP185(:,2)); ylim([2.35 2.65]);xlim([0 24]); ax=gca; ax.XTick=[(0:4:24)];
grid; ylabel('Beta');xlabel('LST (h)');
legend('P0','P5','P180','P185','Location','southeast')
title('Without CMB Removed in the Fit and CMB Not Considered in Beam Adjustment');
text(1,2.4,{['Min F = ' min_f ' MHz'], [Parameters ' Params'], ['Day = ' day]})

%%
% P0 for the 4 cases above
figure(5);clf(5);plot(LST,BetaCMBfitP0(:,2))
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.56 .35 .30 .35])
hold on
plot(LST,BetaNoCMBfitP0(:,2))
plot(LST,BetaCMBfitNoCMBBeamAdjustP0(:,2))
plot(LST,BetaNoCMBfitNoCMBBeamAdjustP0(:,2)); ylim([2.35 2.65]);xlim([0 24]); ax=gca; ax.XTick=[(0:4:24)];
grid; ylabel('Beta');xlabel('LST (h)');
legend('CMB Fit','No CMB Fit','CMB Fit & No CMB Consideration in Beam Adjustment','No CMB Fit & No CMB Consideration in Beam Adjustment','Location','southwest')
title('Beta with P0 various Fit and Beam Scenarios');
text(15,2.6,{['Min F = ' min_f ' MHz'], [Parameters ' Params'], ['Day = ' day]})

set(gcf, 'PaperSize', [7.5 6]); %Set the paper to have width 7.5 and height 6.
set(gcf, 'PaperPosition', [0 0 7.5 6]); %Position plot at left hand corner with width 7.5 and height 6.
savefiledir='C:\Users\Tom\OneDrive\ASU_Thesis\Edges\LOCO Memos\Spectral Index Reports 2015-11 and 2016-06\Figures\';
saveas(gcf, [savefiledir 'CMBVariations_P0_Fmin' min_f 'MHz_Day' day '_' Parameters '_Params.jpg'], 'jpg');


%%
% P185 for the 4 cases above
figure(55);clf(55);plot(LST,BetaCMBfitP185(:,2))
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.70 .35 .30 .35])
hold on
plot(LST,BetaNoCMBfitP185(:,2))
plot(LST,BetaCMBfitNoCMBBeamAdjustP185(:,2))
plot(LST,BetaNoCMBfitNoCMBBeamAdjustP185(:,2)); ylim([2.35 2.65]);xlim([0 24]); ax=gca; ax.XTick=[(0:4:24)];
grid; ylabel('Beta');xlabel('LST (h)');
legend('CMB Fit','No CMB Fit','CMB Fit & No CMB Consideration in Beam Adjustment','No CMB Fit & No CMB Consideration in Beam Adjustment','Location','southwest')
title('Beta with P185 various Fit and Beam Scenarios');
text(15,2.6,{['Min F = ' min_f ' MHz'], [Parameters ' Params'], ['Day = ' day]})

%%
% P0 difference for the 4 cases above
figure(6);clf(6); plot(LST,BetaCMBfitP0(:,2) - BetaNoCMBfitP0(:,2))
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.56 .02 .30 .35])
hold on
plot(LST,BetaCMBfitP0(:,2) - BetaCMBfitNoCMBBeamAdjustP0(:,2))
plot(LST,BetaCMBfitP0(:,2) - BetaNoCMBfitNoCMBBeamAdjustP0(:,2)); ylim([-0.01 +0.03]); xlim([0 24]); ax=gca; ax.XTick=[(0:4:24)];
grid; ylabel('Beta');xlabel('LST (h)');
legend('CMB Fit - No CMB Fit','CMB Fit - CMB Fit but no CMB Consideration in Beam Adjustment','CMB Fit - No CMB Fit & No CMB Consideration in Beam Adjustment','Location','southwest')
title('Beta Differences for P0 Using Various Fit and Beam Scenarios');
text(15,0.025,{['Min F = ' min_f ' MHz'], [Parameters ' Params'], ['Day = ' day]})

set(gcf, 'PaperSize', [7.5 6]); %Set the paper to have width 7.5 and height 6.
set(gcf, 'PaperPosition', [0 0 7.5 6]); %Position plot at left hand corner with width 7.5 and height 6.
savefiledir='C:\Users\Tom\OneDrive\ASU_Thesis\Edges\LOCO Memos\Spectral Index Reports 2015-11 and 2016-06\Figures\';
saveas(gcf, [savefiledir 'CMBDiffs_P0_Fmin' min_f 'MHz_Day' day '_' Parameters '_Params.jpg'], 'jpg');

%%
% P185 difference for the 4 cases above
figure(66);clf(66); plot(LST,BetaCMBfitP185(:,2) - BetaNoCMBfitP185(:,2))
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.70 .02 .30 .35])
hold on
plot(LST,BetaCMBfitP185(:,2) - BetaCMBfitNoCMBBeamAdjustP185(:,2))
plot(LST,BetaCMBfitP185(:,2) - BetaNoCMBfitNoCMBBeamAdjustP185(:,2)); ylim([-0.01 +0.03]); xlim([0 24]); ax=gca; ax.XTick=[(0:4:24)];
grid; ylabel('Beta');xlabel('LST (h)');
legend('CMB Fit - No CMB Fit','CMB Fit - CMB Fit but no CMB Consideration in Beam Adjustment','CMB Fit - No CMB Fit & No CMB Consideration in Beam Adjustment','Location','southwest')
title('Beta Differences for P185 Using Various Fit and Beam Scenarios');
text(15,0.025,{['Min F = ' min_f ' MHz'], [Parameters ' Params'], ['Day = ' day]})

%%
% Case 1 Difference
% Beta with CMB fit and with CMB adjusted beam correction

figure(11);clf(11);plot(LST,BetaCMBfitP0(:,2)-BetaCMBfitP5(:,2));
set(gcf,'Units','normalized','MenuBar','none','ToolBar','none','OuterPosition',[.01 .65 .30 .35])
hold on
plot(LST,BetaCMBfitP0(:,2)-BetaCMBfitP180(:,2))
plot(LST,BetaCMBfitP0(:,2)-BetaCMBfitP185(:,2));
plot(LST,BetaCMBfitP5(:,2)-BetaCMBfitP185(:,2));ylim([-0.004 0.008]); xlim([0 24]); ax=gca; ax.XTick=[(0:4:24)];
grid; ylabel('Beta');xlabel('LST (h)');
legend('P0 - P5','P0 - P180','P0 - P185','P5 - P185','Location','southeast');
title('Difference Between Polarization Angles Full CMB Considerations');
text(15,7e-3,{['Min F = ' min_f ' MHz'], [Parameters ' Params'], ['Day = ' day]})

set(gcf, 'PaperSize', [7.5 6]); %Set the paper to have width 7.5 and height 6.
set(gcf, 'PaperPosition', [0 0 7.5 6]); %Position plot at left hand corner with width 7.5 and height 6.
savefiledir='C:\Users\Tom\OneDrive\ASU_Thesis\Edges\LOCO Memos\Spectral Index Reports 2015-11 and 2016-06\Figures\';
saveas(gcf, [savefiledir 'FullCMB_PDiffs_Fmin' min_f 'MHz_Day' day '_' Parameters '_Params.jpg'], 'jpg');

%%
figure(1); figure(2);figure(5); % regain focus/overlap on screen

