close all; clear all;
% global pets 

load('results_Danio_rerio_limited.mat')
[data, auxData, metaData, txtData, weights] = mydata_Danio_rerio_limited; 
% [par, metaPar, txtPar] = pars_init_Danio_rerio_limited(metaData);
par
[prdData_limited, info] = predict_Danio_rerio_limited(par,data,auxData);

load('results_Danio_rerio_moderate.mat')
[data, auxData, metaData, txtData, weights] = mydata_Danio_rerio_moderate; 
% [par, metaPar, txtPar] = pars_init_Danio_rerio_moderate(metaData);
[prdData_moderate, info] = predict_Danio_rerio_moderate(par,data,auxData);

% load('results_Danio_rerio_full.mat')
% [data, auxData, metaData, txtData, weights] = mydata_Danio_rerio_full; 
% % [par, metaPar, txtPar] = pars_init_Danio_rerio_full(metaData);
% [prdData_full, info] = predict_Danio_rerio_full(par,data,auxData);

[prdData_limited.Ri_28;prdData_moderate.Ri_28]

figure(1)
x = data.tL_BarrFern2010(:,1);
plot(x, prdData_limited.tL_BarrFern2010,'r', 'LineWidth',3.5), hold on
plot(x, prdData_moderate.tL_BarrFern2010,'g', 'LineWidth',3.5), hold on
% plot(x, prdData_full.tL_BarrFern2010,'b', 'LineWidth',3.5), hold on
plot(x, data.tL_BarrFern2010(:,2),'*', 'LineWidth',3), hold on
% legend('low', 'moderate', 'high', 'Location', 'best', 'FontSize',20)
legend('low', 'moderate', 'Location', 'best', 'FontSize',20)
xlabel('time, d', 'FontSize',18)
ylabel('length, mm', 'FontSize',18)

figure(2)
x = data.tW_BarrFern2010(:,1);
plot(x, prdData_limited.tW_BarrFern2010,'r', 'LineWidth',3.5), hold on
plot(x, prdData_moderate.tW_BarrFern2010,'g', 'LineWidth',3.5), hold on
% plot(x, prdData_full.tW_BarrFern2010,'b', 'LineWidth',3.5), hold on
plot(x, data.tW_BarrFern2010(:,2),'*', 'LineWidth',3), hold on
% legend('low', 'moderate', 'high', 'Location', 'best', 'FontSize',20)
legend('low', 'moderate', 'high', 'Location', 'best', 'FontSize',20)
xlabel('time, d', 'FontSize',18)
ylabel('wet weight, mg', 'FontSize',18)

figure(3)
x = data.tJO(:,1);
plot(x, prdData_limited.tJO,'r', 'LineWidth',3.5), hold on
plot(x, prdData_moderate.tJO,'g', 'LineWidth',3.5), hold on
% plot(x, prdData_full.tJO,'b', 'LineWidth',3.5), hold on
plot(x, data.tJO(:,2),'*', 'LineWidth',3), hold on
% legend('low', 'moderate', 'high', 'Location', 'best', 'FontSize',20)
legend('low', 'moderate', 'high', 'Location', 'best', 'FontSize',20)
xlabel('time, d', 'FontSize',18)
ylabel('oxygen consumption, \mumol/g/h', 'FontSize',18)
