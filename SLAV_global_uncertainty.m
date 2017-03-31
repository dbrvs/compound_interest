function [RHOh,RHOp]=SLAV_uncertainty(numpatient,pot,tf)
% integrate the stochastic SLAV differential equations with LHS sampling
% DBR 4/16

%numpatient = 100; % number of patients


%variables to allow to be variable, some need to over multi logs

%      eps   R0 fn  aLc aLe aLn  g  dS  dA  
lb = [ 0.7   0  1   1   3.8 .15  3  .02 0.8 ]; % lower bound
ub = [ 0.99  1  10  2   5.7 .23  23 .4  1.2 ]; % upper bound

%        tau  bt0  xi  aS  p   L0 thL 
lblog = [ -5  -5   -7  1   2   4  -4 ]; % lower bound for multilog params
ublog = [ -3  -3   -3  3   4   8  -3 ]; % upper bound for multilog

param1 = length(lb); % number of parameters
xn1 = lhsdesign(numpatient,param1); % generate normalized design
lm1 = bsxfun(@plus,lb,bsxfun(@times,xn1,(ub-lb)));

param2 = length(lblog); % number of parameters
xn2 = lhsdesign(numpatient,param2); % generate normalized design
lm2 = 10.^bsxfun(@plus,lblog,bsxfun(@times,xn2,(ublog-lblog)));

lm = [lm1 lm2 1+lhsdesign(numpatient,1)*(pot-1) ];

% now the matrix of parameters is complete, calculate a few other things
% 1   2  3   4   5   6   7  8  9   10  11  12 13 14 15 16
% eps R0 fn  aLc aLe aLn g  dS dA  tau bt0 xi aS p  L0 thL 

n   = lm(:,14)./lm(:,7);
S0  = lm(:,13)./lm(:,8);
Lf  = 1-lm(:,10).*(1+lm(:,12)./lm(:,16));
eps = lm(:,1);
bte = lm(:,11).*(1-lm(:,1));
R0  = bte.*S0.*n./lm(:,9).*Lf;
L0  = lm(:,15)/5e6;

thL = -lm(:,16); 
xi  = lm(:,12);
tau = lm(:,10);

fT  = [1-lm(:,3)/100-0.4 repmat(0.4,numpatient,1) lm(:,3)/100]; %fractions, assume Tem~0.4 always

aL  = [lm(:,4) lm(:,5) lm(:,6)]/100;

dL  = bsxfun(@minus,aL,thL+xi);

potency=lm(:,17); %compound interest potency

%% solve the odes for "good patients" with R0<1

%tf = 100;  %years to simulate
goodpatients = find(R0<1);

leg={};

for j=1:length(goodpatients)
    
in=goodpatients(j);

ops = odeset('RelTol', 1e-10); %tolerances for ODE solver

X0=[S0(in) L0(in)*fT(in,1) L0(in)*fT(in,2) L0(in)*fT(in,3) 3e-3];      

%solve the ODE twice, once for upper and lower range
ts=linspace(0,tf,1e3)*365;

aL_CI=aL(in,:)/potency(in,:); %compound interest cure potency

%aS,dS,bt,tau,aL,dL,xi,fT,dA,n

[t, Y]=ode23s(@subsetSLAV,ts,X0,ops,...
        lm(in,13),lm(in,8),bte(in),lm(in,10),aL_CI,dL(in,:),xi(in),fT(in,:),lm(in,9),n(in));

ty(:,j)=t/365; 

Ltot(:,j)=(Y(:,2)+Y(:,3)+Y(:,4))*5e6; %total number of latent cells

V(:,j) = Y(:,5)*n(in)*1e3;

leg{j}=num2str(R0(in));
% 
 index2pink=find(Ltot(:,j)<1e4);
 index2hill=find(Ltot(:,j)<1e2);
 if isempty(index2pink)~=0 || isempty(index2hill)~=0
     hillt(j)=NaN;
     pinkt(j)=NaN;
 else
     hillt(j)=ty(index2hill(1));
     pinkt(j)=ty(index2pink(1));
 end

end

figure(1)
 semilogy(ty,Ltot)
 xlim([0,50])
 %ylim([100,1e8])
 figuresize(4,6,'inches')
 children = findall(gcf); for k = 1:numel(children); try set(children(k), 'FontSize',14); end; end
 print('SLAV_uncert','-dpdf')

figure(2)
figuresize(6,6,'inches')
varz=[tau eps R0 L0*5e6 -thL potency fT(:,3)];
varzname={'\tau','\epsilon^{ART}', 'R_0^{ART}','L_0 (cells)','-\theta_L (1/day)','\epsilon^{AP}','L_n(0)/L_0'};
for i=1:size(varz,2)
    [RHOh(i),~] = corr(varz(goodpatients,i),hillt');
    [RHOp(i),~] = corr(varz(goodpatients,i),pinkt');
    %plot scatterplots of data correlations
    subplot(4,2,i)
    scatter(varz(goodpatients,i),pinkt,10,'black','filled')
    hold on
    scatter(varz(goodpatients,i),hillt,10,'red','filled')
    set(gca,'yscale','log')
    hold off
    xlim([min(varz(goodpatients,i)) max(varz(goodpatients,i))])
    ylim([0.1 100])
        set(gca,'ytick',[1 10 100]);
        set(gca,'yticklabel',['10^0'; '10^1'; '10^2']);
    xlabel(varzname{i})
%    title(['\rho=' num2str(RHOp,2) ' (p=' num2str(PVAL,2) ')']);

end

children = findall(gcf); for k = 1:numel(children); try set(children(k), 'FontSize',10); end; end
print('SLAV_corr','-dpdf')

%subplot(339)
figure(3)
figuresize(4,6,'inches')
%ll=[RHOp;RHOh];
[rhops,indx]=sort(RHOp);

ll=[rhops;RHOh(indx)];

b=barh(1:length(ll),ll');
b(2).FaceColor = 'red';
b(1).FaceColor = 'black';
xlim([-1,1])
ylim([0,length(ll)+1])
%set(gca,'yticklabel','');
set(gca,'yticklabel',varzname(indx));
legend('P1','Hc')
legend('Location','SouthEast')
%set(gca,'yticklabel','');

children = findall(gcf); for k = 1:numel(children); try set(children(k), 'FontSize',14); end; end
print('SLAV_corr_bar','-dpdf')

%close all

