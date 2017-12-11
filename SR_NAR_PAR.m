close all;clear all;
tsim=100;
%% Simple Regulation Model
SR = sbiomodel('SimpleRegulation');
Dil = 0.25;     % Deradation/dilution rate 
ActY = 0.1;     % Activation strength of Y
Tr = 0.5;       % Translation rate

X=0;Y=0.5;XmRNA = 0;
speciesObj = addspecies(SR,'X',X);
speciesObj = addspecies(SR,'Y',Y);
set (speciesObj, 'ConstantAmount', true);

speciesObj = addspecies(SR,'XmRNA',XmRNA);

% X degradation/dilution
V2 = 'K2*X';
r = addreaction(SR, 'X -> null');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K2', Dil);
set (r, 'ReactionRate', V2);

% Translation to X
V3 = 'XmRNA*K3'; 
r = addreaction(SR, 'XmRNA -> X');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K3', Tr);
set (r, 'ReactionRate', V3);

% Y Induces XmRNA transcription
V5 = 'Y*K5'; 
r = addreaction(SR, 'Y -> XmRNA');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K5', ActY);
set (r, 'ReactionRate', V5);

% XmRNA Dgradation
V6 = 'XmRNA*K6';
r = addreaction(SR, 'XmRNA -> null');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K6', Dil);
set (r, 'ReactionRate', V6);

% Y Degradation
V7 = 'K2*Y';
r = addreaction(SR, 'Y -> null');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K2', Dil);
set (r, 'ReactionRate', V7);

csObj = getconfigset(SR,'active'); set(csObj,'Stoptime',tsim);
[t1,x1,names1] = sbiosimulate(SR);


%% Negative Autoregulation Model
NAR = sbiomodel('NegativeAutoregulation');

Rep = 0.5;  % Repression strength
Pro = 1;    % Promoter strength

X=0; XmRNA = 0;
speciesObj = addspecies(NAR,'X',X);
speciesObj = addspecies(NAR,'XmRNA',XmRNA);

% Translation of X
V1 = 'XmRNA*K1';
r = addreaction(NAR, 'XmRNA -> X');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K1', Tr);
set (r, 'ReactionRate', V1);

% X degradation
V2 = 'K2*X';
r = addreaction(NAR, 'X -> null');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K2', Dil);
set (r, 'ReactionRate', V2);

% XmRNA degradation
V3 = 'XmRNA*K3'; 
r = addreaction(NAR, 'XmRNA -> null');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K3', Dil);
set (r, 'ReactionRate', V3);

% X negatively acting on its own transcription
V4 = 'K4-X*K5'; 
r = addreaction(NAR, 'null -> XmRNA');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K4', Pro);
p = addparameter(k, 'K5', Rep);
set (r, 'ReactionRate', V4);


csObj = getconfigset(NAR,'active'); set(csObj,'Stoptime',tsim);
[t2,x2,names2] = sbiosimulate(NAR);

% Testing different repression parameters
f2= figure;
Rep = [0.01,0.1,1];
for i = 1:3
    NAR.Reactions(4).KineticLaw.Parameters(2).Value = Rep(i);
    [t2a,x2a,names2a] = sbiosimulate(NAR);
    plot(t2a,x2a(:,1)/x2a(end,1),'LineWidth',2,'LineSmoothing','on');
    hold on;
end
set(gca,'LineWidth',2); set(gca,'fontsize',16);
h=legend('0.01','0.1','1','Location','southeast');
v = get(h,'title'); set(v,'string','Repression parameter');
ylabel('X/Xst');xlabel('time');xlim([0 50]);
print('f2','-dpng','-r300')

%% Positive Autoregulation Model
PAR = sbiomodel('PositiveAutoregulation');
X=0; XmRNA = 0;
speciesObj = addspecies(PAR,'X',X);
speciesObj = addspecies(PAR,'XmRNA',XmRNA);

WeakPro = 1;
Act = 0.2;

V1 = 'XmRNA*K1'; 
r = addreaction(PAR, 'XmRNA -> X');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K1', Tr);
set (r, 'ReactionRate', V1);

V2 = 'K2*X';
r = addreaction(PAR, 'X -> null');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K2', Dil);
set (r, 'ReactionRate', V2);

V3 = 'XmRNA*K3';
r = addreaction(PAR, 'XmRNA -> null');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K3', Dil);
set (r, 'ReactionRate', V3);

V4 = 'X*K4+K5'; 
r = addreaction(PAR, 'null -> XmRNA');
k = addkineticlaw(r, 'Unknown');
p = addparameter(k, 'K4', Act);
p = addparameter(k, 'K5', WeakPro);
set (r, 'ReactionRate', V4);

csObj = getconfigset(PAR,'active'); set(csObj,'Stoptime',tsim+100);
[t3,x3,names3] = sbiosimulate(PAR);
Act = [0.01,0.1,0.3];
f3 = figure;
for i = 1:3
    % Only Act change the timing before reaching steadystate
    PAR.Reactions(4).KineticLaw.Parameters(1).Value = Act(i);
    [t3a,x3a,names3a] = sbiosimulate(PAR);
    plot(t3a,x3a(:,1)/x3a(end,1),'LineWidth',2,'LineSmoothing','on');
    hold on;    
end
set(gca,'LineWidth',2); set(gca,'fontsize',16);
ylabel('X/Xst');xlabel('time');
h = legend('0.01','0.1','0.3','Location','southeast');
v = get(h,'title'); set(v,'string','Activation parameter');
xlim([0 50])
print('f3','-dpng','-r300')

%% Plot for the 3 models comparison
X1 = x1(:,1)/x1(end,1); X2 = x2(:,1)/x2(end,1); X3 = x3(:,1)/x3(end,1);
f5 = figure;
plot(t1,X1(:,1),'b','LineWidth',2,'LineSmoothing','on'); 
set(gca,'LineWidth',2); set(gca,'fontsize',16); hold on;
plot(t2,X2(:,1),'g','LineWidth',2,'LineSmoothing','on'); 
set(gca,'LineWidth',2); set(gca,'fontsize',16); hold on;
plot(t3,X3(:,1),'r','LineWidth',2,'LineSmoothing','on'); 
set(gca,'LineWidth',2); set(gca,'fontsize',16); hold on;
legend('Simple regulation','Negative autoregulation',...
    'Positive autoregulation','Location','southeast')

xlim([0 50])
ylabel('X/Xst');xlabel('time');print('f5','-dpng','-r300');