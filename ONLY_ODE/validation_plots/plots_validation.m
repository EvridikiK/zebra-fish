close all; clear all;
% global pets 

load('results_Danio_rerio_limited.mat')
[data, auxData, metaData, txtData, weights] = mydata_Danio_rerio_limited; 
[par, metaPar, txtPar] = pars_init_Danio_rerio_limited(metaData);
% par
[prdData_limited, info] = predict_Danio_rerio_limited(par,data,auxData);

load('results_Danio_rerio_moderate.mat')
[data, auxData, metaData, txtData, weights] = mydata_Danio_rerio_moderate; 
[par, metaPar, txtPar] = pars_init_Danio_rerio_moderate(metaData);
[prdData_moderate, info] = predict_Danio_rerio_moderate(par,data,auxData);

load('results_Danio_rerio_full.mat')
[data, auxData, metaData, txtData, weights] = mydata_Danio_rerio_full; 
[par, metaPar, txtPar] = pars_init_Danio_rerio_full(metaData);
[prdData_full, info] = predict_Danio_rerio_full(par,data,auxData);

[prdData_limited.Ri_28;prdData_moderate.Ri_28;prdData_full.Ri_28]

figure(1)
x = data.tL_BarrFern2010(:,1);
plot(x, prdData_limited.tL_BarrFern2010, 'Color', [0, 0.4470, 0.7410], 'LineWidth',3.5), hold on
plot(x, prdData_moderate.tL_BarrFern2010, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth',3.5), hold on
plot(x, prdData_full.tL_BarrFern2010,'Color', [0.9290, 0.6940, 0.1250], 'LineWidth',3.5), hold on
plot(x, data.tL_BarrFern2010(:,2),'*','Color', [0.4660, 0.6740, 0.1880], 'LineWidth',3), hold on
legend('low', 'moderate', 'high', 'Location', 'best', 'FontSize',20,'Box','off')
xlabel('time, d', 'FontSize',18)
ylabel('length, mm', 'FontSize',18)
ax = gca;
ax.FontSize = 20;

figure(2)
x = data.tW_BarrFern2010(:,1);
plot(x, prdData_limited.tWw_BarrFern2010, 'Color', [0, 0.4470, 0.7410], 'LineWidth',3.5), hold on
plot(x, prdData_moderate.tWw_BarrFern2010, 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth',3.5), hold on
plot(x, prdData_full.tWw_BarrFern2010, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth',3.5), hold on
plot(x, data.tW_BarrFern2010(:,2), '*','Color', [0.4660, 0.6740, 0.1880], 'LineWidth',3), hold on
legend('low', 'moderate', 'high', 'Location', 'best', 'FontSize',20,'Box','off')
xlabel('time, d', 'FontSize',18)
ylabel('wet weight, mg', 'FontSize',18)
ax = gca;
ax.FontSize = 20;

figure(3)
x = data.tJO(:,1);
plot(x, prdData_limited.tJO, 'Color', [0, 0.4470, 0.7410], 'LineWidth',3.5), hold on
plot(x, prdData_moderate.tJO,'Color', [0.8500, 0.3250, 0.0980], 'LineWidth',3.5), hold on
plot(x, prdData_full.tJO,'Color', [0.9290, 0.6940, 0.1250], 'LineWidth',3.5), hold on
plot(x, data.tJO(:,2), '*','Color', [0.4660, 0.6740, 0.1880], 'LineWidth',3), hold on
legend('low', 'moderate', 'high', 'Location', 'best', 'FontSize',20,'Box','off')
xlabel('time, d', 'FontSize',18)
ylabel('oxygen consumption, \mumol/g/h', 'FontSize',18)
ax = gca;
ax.FontSize = 20;
