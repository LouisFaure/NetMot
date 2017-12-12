load('C1FFL_I1FFL');
%%
for i = 1:2
    if i==1     % CFFL1
        C1FFL_I1FFL.Reactions(2).Active = 1;
        C1FFL_I1FFL.Reactions(3).Active = 0;
        idx = 'a';
    elseif i==2 % IFFL1 
        C1FFL_I1FFL.Reactions(2).Active = 0;
        C1FFL_I1FFL.Reactions(3).Active = 1;
        idx = 'b';
    end
    csObj = getconfigset(C1FFL_I1FFL,'active'); set(csObj,'Stoptime',20);
    eval(['[t',idx,',x',idx,',names] = sbiosimulate(C1FFL_I1FFL,csObj);'])
    eval(['t = t',idx,';']); eval(['x = x',idx,';']);
    
    % Produce plot of Sx induction
    if i==1
        f4 = figure;
        plot(t,x(:,6),'k','LineWidth',2,'LineSmoothing','on');
        ylabel('Sx');
        set(gca,'LineWidth',2); set(gca,'fontsize',16);
        legend(names(6)); ylim([0 1.5]);f4.Position(3) = f4.Position(3)*2;
        print('f4','-dpng','-r300');
    end
    
    f4_2 = figure;
    plot(t,ones(size(t)),'--k','LineWidth',2,'LineSmoothing','on'); hold on;
    plot(t,x(:,[4,5])./x(find(x(:,2)==1, 1, 'last' ),[4,5])...
        ,'LineWidth',2,'LineSmoothing','on'); 

    set(gca,'LineWidth',2); set(gca,'fontsize',16); 
    legend('Steady state',names{3},names{5}); ylim([0 2]);
    xlabel('time'); ylabel('SP/SPst'); f4_2.Position(3) = f4_2.Position(3)*2;
    print(['f4',idx],'-dpng','-r300');
end


%% Differentiation
clear all; close all;
load('Diff');
csObj = getconfigset(Diff,'active'); set(csObj,'Stoptime',20);
[t,x,names] = sbiosimulate(Diff,csObj);


%% Plotting induction and Z pulses
f6 = figure;
subplot(2,1,1);
plot(t,x(:,6),'k','LineWidth',2,'LineSmoothing',   'on'); 
set(gca,'LineWidth',2); set(gca,'fontsize',16); hold on;
legend(names(6));ylim([0 1.5]);ylabel('amount');xlim([0 15]);
subplot(2,1,2);
plot(t,x(:,[5,10,11]),'LineWidth',2,'LineSmoothing',   'on'); 
set(gca,'LineWidth',2); set(gca,'fontsize',16); hold on;
legend(names([5,10,11]));ylabel('amount');xlabel('time'); xlim([0 15]);
f6.Position = [1 1 f6.Position(3:4).*2];
print('f5','-dpng','-r300');