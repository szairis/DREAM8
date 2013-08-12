% PRODUCE VISUAL OUTPUT FOR THE PRESCRIBED EXPERIMENT.
% Inputs: refer to the curves.m function.
% © Chris Oates 2011.

function Investigate(Generator,its_max,t,snr_cell,snr_meas,mod_S)

% obtain data
[data,genes,TPR,FPR,AUR,P,R,schemes,table] = curves(Generator,its_max,t,snr_cell,snr_meas,mod_S);

% plot typical simulation data
figure
datasets = size(data,3);
grid = ceil(sqrt(datasets));
for i = 1:datasets
    subplot(grid,grid,i)
    scale = 1.2*max(max(data(:,:,1)));
    set(0,'DefaultAxesColorOrder',[0 0 0]);
    c = {'s','d','^','v','>','<','+','o','*','x','p','hexagram'};
    set(0,'DefaultAxesLineStyleOrder',c);
    set(gca,'FontSize',24/grid)
    plot(t,data(:,1:length(t),i)','MarkerSize',15/grid,'LineWidth',3/grid)
    hold on
    plot(t,data(:,1:length(t),i)','--','LineWidth',3/grid)
    hold off
    legend(genes,'FontSize',24/grid,'Interpreter','none')
    xlim([t(1) t(end)])
    ylim([0 scale])
    xlabel('Time (mins)','FontSize',24/grid)
    ylabel('Expression Level (arbitrary units)','FontSize',24/grid)
end
drawnow

% plot ROC curves
figure
D = length(schemes);
set(0,'DefaultAxesColorOrder',colormap(varycolor(D)))
set(gca,'FontSize',24)
plot(FPR',TPR','--','LineWidth',3)
set(gca,'XTick',[0 1])
set(gca,'YTick',[0 1])
xlim([0 1])
ylim([0 1])
xlabel('False Positive Rate','FontSize',24)
ylabel('True Positive Rate','FontSize',24)
legend(schemes,'FontSize',9,'Interpreter','none')
drawnow

% plot PR curves
figure
set(gca,'FontSize',24)
plot(R',P','--','LineWidth',3)
set(gca,'XTick',[0 1])
set(gca,'YTick',[0 1])
xlim([0 1])
ylim([0 1])
xlabel('Recall','FontSize',24)
ylabel('Precision','FontSize',24)
legend(schemes,'FontSize',9,'Interpreter','none')
drawnow

% plot AUR scores
figure('Position',[100 50 600 650])
h1 = subplot(1,2,1);
p1 = [0.08 0.08 0.42 0.9];
set(h1,'pos',p1);
b = boxplot(AUR,'orientation','horizontal');
set(b(:,1:end),'LineWidth',3);
set(findobj(gca,'Type','text'),'FontSize',24)
set(gca,'YTick',[0.5:1:31.5])
set(gca,'YTickLabel',{' '})
set(gca,'YGrid','on','LineWidth',3)
ylabel('Area Under ROC Curve','FontSize',24)
xlim([0 1])
set(gca,'FontSize',24)
set(gca,'XTick',[0 0.5 1])
hold on
plot([0.5 0.5],ylim,'--','LineWidth',3)
hold off
h2 = subplot(1,2,2);
p2 = [0.5 0.08 0.42 0.9];
set(h2,'pos',p2);
clims = [-5 2];
table = table(end:-1:1,:);
imagesc(table,clims)
colormap(gray)
[I,J] = size(table);
textHandles = zeros([I,J])';
for i = 1:I
    for j = 1:J 
        if j == 1
            if table(i,j) == 0
                e = 'AIC';
            else
                e = 'Bayes';
            end
        elseif j == 2
            if table(i,j) == 0
                e = 'Std';
            else
                e = 'Quad';
            end
        elseif j == 3
            if table(i,j) == 0
                e = 'None';
            else
                e = 'Used';
            end
        elseif j == 4
            if table(i,j) == -1
                e = 'x';
            elseif table(i,j) == 0;
                e = '0';
            elseif table(i,j) == 1;
                e = '1';
            else
                e = '2';
            end
        end
        textHandles(j,i) = text(j,i,e,'horizontalAlignment','center','fontsize',10,'clipping','on','visible','off');
    end
end
set(textHandles,'visible','on')
set(gca,'YGrid','on','LineWidth',3)
set(gca,'YTickLabel',{' '})
set(gca,'XTick',[1 2 3 4])
set(gca,'FontSize',14)
set(gca,'XTickLabel',{'Inf.','D. Mat','Lag','Var. Fun.'})
set(gca,'YTick',[0.5:1:31.5])

end
