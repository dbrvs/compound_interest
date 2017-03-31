%%% code that makes a bunch of plots based on the decoupled compound
%%% interest formula
%%% dbr 10/15/15

%% Figure 2a -- lines comparison

clear all

aL  = 0.015;  %growth rate latent cells [1/day]
xi  = 5.7e-5; %reactivation rate in Luo patient specific 2e-6->1e-3
thL = -5.2e-4; %net clearance rate

dL = aL-thL-xi; 

thL=aL-dL-xi;
tt=linspace(0,100,1e3); %timeseries in years
t=tt.*365; %convert to days

figuresize(3,3,'inches')

Hc = ones(length(t),1)*100; %Hill cure in half threshold
H1 = ones(length(t),1)*500; %Hill 1 year rebound threshold
D1 = ones(length(t),1)*1.5e4; %davenport 1 year rebound threshold

% one time therapy affects L0

L= 1e6*exp(thL*t); 
%subplot(121)
semilogy(tt,L,tt,L/10,'--',tt,L/100,':',tt,L/1e3,'-.','Linewidth',3)
l=legend('Natural clearance','10-fold decrease L(t=0)','100-fold','1,000-fold');
set(l,'location','SouthWest');
set(l,'FontSize',10);

xlim([0,30]); 
set(gca,'XTick',0:5:30)

ylim([1,2e6])
set(gca,'YTick',logspace(0,7,8))

hold on
semilogy(tt,D1,tt,H1,tt,Hc,'color','black','linestyle','-','linewidth',2)
hold off

%cs = findall(gcf); for k = 1:numel(cs); try set(cs(k), 'FontSize',14); end; end
print('PIC_onetime','-dpdf','-r600')


%% Figure 2b -- lines comparison

% continuous time therapy affects theta_L
figuresize(5,5,'inches')


thL_c = [thL thL*1.5 thL*3 thL*10];
L_c = 1e6*exp([thL_c(1)*t; thL_c(2)*t; thL_c(3)*t; thL_c(4)*t]);
%subplot(122)
semilogy(tt,L_c(1,:),tt,L_c(2,:),'--',tt,L_c(3,:),':',tt,L_c(4,:),'-.','Linewidth',3)
legend('Natural clearance','1.5-fold increase \theta_L','3-fold','10-fold')
legend('location','SouthWest')
xlim([0.05,30])
ylim([1,2e6])
%set(gca,'YTickLabel','')
hold on
semilogy(tt,D1,tt,H1,tt,Hc,'color','black','linestyle','-','linewidth',2)
hold off

cs = findall(gcf); for k = 1:numel(cs); try set(cs(k), 'FontSize',14); end; end

xlim([0,30]); 
set(gca,'XTick',0:5:30)

ylim([1,2e6])
set(gca,'YTick',logspace(0,7,8))

print('PIC_continuous','-dpdf','-r600')

%close all

%% Figure 2c heat maps duration vs potency

clear all

aL  = 0.015;  %growth rate latent cells [1/day]
xi  = 5.7e-5; %reactivation rate in Luo patient specific 2e-6->1e-3
thL = -5.2e-4; %net clearance rate
dL  = aL-thL-xi; 

% changing xi for kick and kill

for i=1:156
for tf=1:156 %[weeks]
    XI(i)=xi*i;
    thL(i)=aL-dL-XI(i);
    numleft(tf,i) = 1e6*exp(tf*7*thL(i));
    pctleft(tf,i) = exp(tf*7*thL(i));
    tw(tf)=tf;
end
end

figuresize(3,3,'inches')
%subplot(121)
[hC1 hC1] = contourf(1:i,tw,log10(numleft),20);
set(gca,'XScale','log');
caxis([1,6])
set(hC1,'LineStyle','none');
hold on
contour(1:156,tw,log10(numleft),[2 2],'color','k','linewidth',3)
contour(1:156,tw,log10(numleft),[2.7 2.7],'color','k','linewidth',3)
contour(1:156,tw,log10(numleft),[4.1 4.1],'color','k','linewidth',3)
hold off

% inset the contour lines

%axes('Position',[.17 .2 .1 .2])
%box on
%contourf(1:10,1:10,log10(numleft(1:10,1:10)))
% = contourf(peaks);
%caxis([1,6])

colormap('jet')     
cs=findall(gcf); for k = 1:numel(cs); try set(cs(k), 'FontSize',14); end; end

xlim([0,150]); 
set(gca,'YTick',0:25:150)

ylim([0,150])
set(gca,'XTick',logspace(0,2,3))

print('PIC_mapLRA','-dtiff','-r600')

%% Figure 2d heat maps duration vs potency
% antirpoliferative therapy

for i=1:156
for tf=1:156 %[weeks]
    A(i)=aL/i;
    thL(i)=A(i)-dL-xi;
    numleft(tf,i) = 1e6*exp(tf*7*thL(i));
    tw(tf)=tf;
end
end

%subplot(122)

figuresize(3,3,'inches')
[hC2 hC2] = contourf(1:156,tw,log10(numleft),20);
set(gca,'XScale','log'); 
set(hC2,'LineStyle','none');
hold on
contour(1:156,tw,log10(numleft),[2 2],'color','k','linewidth',3)
contour(1:156,tw,log10(numleft),[2.7 2.7],'color','k','linewidth',3)
contour(1:156,tw,log10(numleft),[4.1 4.1],'color','k','linewidth',3)
hold off

%p2SizeI = get(gca, 'Position'); %used to adjust for colorbar changing the size!

%colorbar
%set(gca, 'Position', p2SizeI-[0.05 0 0 0]); % Can also use gca instead of h2 if h2 is still active.

cs=findall(gcf); for k = 1:numel(cs); try set(cs(k), 'FontSize',14); end; end

%axis tight
caxis([1,6])
%set(gca,'ZTick',1:5,'ZTickLabel',[1e1 1e2 1e3 1e4 1e5])
colormap('jet')     
%h = colorbar;
%set( h, 'YDir', 'reverse' );

xlim([0,150]); 
set(gca,'YTick',0:25:150)

ylim([0,150])
set(gca,'XTick',logspace(0,2,3))

print('PIC_mapAP','-dtiff','-r600')


%% Fig 2e comparing directly potency of LRA vs AP

figuresize(3,3,'inches')

for i=1:156
for j=1:156 %[weeks]
    A(i)=aL/i;
    XI(j)=xi*j;
    thL=A(i)-dL-XI(j);
    numleft(j,i) = 1e6*exp(7*70*thL);
end
end

[hC2 hC2] = contourf(1:156,1:156,log10(numleft),20);
set(gca,'XScale','log'); 
set(gca,'YScale','log');
set(hC2,'LineStyle','none');
hold on
contour(1:156,1:156,log10(numleft),[2 2],'color','k','linewidth',3)
contour(1:156,1:156,log10(numleft),[2.7 2.7],'color','k','linewidth',3)
%contour(1:156,tw,log10(numleft),[2.6 2.6],'color','k','linewidth',3)
contour(1:156,1:156,log10(numleft),[4.1 4.1],'color','k','linewidth',3)
hold off

%p2SizeI = get(gca, 'Position'); %used to adjust for colorbar changing the size!

%colorbar
%set(gca, 'Position', p2SizeI-[0.05 0 0 0]); % Can also use gca instead of h2 if h2 is still active.

cs=findall(gcf); for k = 1:numel(cs); try set(cs(k), 'FontSize',14); end; end

%axis tight
caxis([1,6])
%set(gca,'ZTick',1:5,'ZTickLabel',[1e1 1e2 1e3 1e4 1e5])
colormap('jet')     
%h = colorbar;
%set( h, 'YDir', 'reverse' );

xlim([0,150]); 
set(gca,'YTick',logspace(0,2,3))

ylim([0,150])
set(gca,'XTick',logspace(0,2,3))

print('PIC_map3','-dtiff','-r600')

%% Figure 4 logistic model of latency


aL   = 0.015;  %growth rate latent cells [1/day]
dL   = 0.0153; %death rate latent cells [1/day]
xi   = 5.7e-5; %reactivation rate in Luo patient specific 2e-6->1e-3

Lmax=1e6;

a=aL-dL-xi;
b=aL/Lmax;

L0=1e6;

c1=L0/(b*L0-a);

t=(0:.1:70)*365;

Ll=a*c1*exp(a*t)./(b*c1*exp(a*t)-1);
L=L0*exp(a*t);

figuresize(4,4,'inches')
semilogy(t/365,L,'--',t/365,Ll,'linewidth',3)
ylabel('latent cells remaining')
xlabel('time (years)')
legend('Typical','Logistic growth')
legend('Location','SouthWest')
legend('boxoff')

cs=findall(gcf); for k = 1:numel(cs); try set(cs(k), 'FontSize',14); end; end


print -dpdf PIC_logisticLAT.pdf


%% jacobian eigenvalues
clear all

%model parameters from nature methods
tau  = 1e-4;    %fraction of infecteds that go latent []
aL   = 0.015;  %growth rate latent cells [1/day]
dL   = 0.015509; %death rate latent cells [1/day]
aA   = 0;      %growth rate active cells [1/day]
dA   = 1;      %death rate active cells [1/day] from Luo total change is -0.18->-2.3 
beta = 1e-4;   %infectivity [uL/virion/day]
xi   = 5.7e-5; %reactivation rate in Luo patient specific 2e-6->1e-3
aS   = 300;    %constant growth rate of T cells 35-760 [cells/uL/day], I converted using body 6L to get [cells/day]
dS   = 0.2;    %susceptible death rate [1/day] Luo big range 0.045->0.45
g    = 23;     %clearance rate Luo 18.8 [1/day]
p    = 1e3;      %viral expansion rate Luo 2.4-9.8 [virions/cell]
thL  = aL-dL-xi;

fL=1-(1+xi/thL)*tau;

cep = 1-dS*g*dA/p/beta/aS/fL;

for i=1:90
    
    ep(i)=i/100;
    
betaE=beta*(1-ep(i));

S=g*dA/betaE/p/fL;
L=tau/thL*(g*dA*dS/betaE/p/fL-aS);
A=aS*fL/dA-g*dS/p/betaE;
V=p*aS*fL/g/dA-dS/betaE;

Jeq =[[-dS-betaE*V 0 0 -betaE*S];[tau*betaE*V thL 0 tau*betaE*S];...
     [(1-tau)*betaE*V xi aA-dA (1-tau)*betaE*S];[0 0 p -g]];
 
eD(i,:)=eig(Jeq);

slav(i,:)=[S L A V];

ctime(i)=1./min(abs(eD(i,:)))/365;

end
%
for i=1:4
subplot(4,1,i)
plot(ep,real(eD(:,i)))
hold on
plot(ep,imag(eD(:,i)))
hold off
end

print('PIC_eigenvals_ep','-dpdf','-r600')

%% figure 3 plots of heterogeneous latent pools

thL=-5.2e-4;
xi = 1.1e-5;

aTn = 0.002;
dTn = aTn-thL-xi;

aTc = 0.015;
dTc = aTc-thL-xi;

aTe = 0.047;
dTe = aTe-thL-xi;

thc = aTc/10-dTc-xi;
the = aTe/10-dTe-xi;
thn = aTn/10-dTn-xi;

%full calculation

Tn=linspace(0,40,100); %percent naive
Te=linspace(0,40,100); %percent em
 
t=linspace(0,10,10000);

tt=t*365; %time in days

L=zeros(length(Tn),length(Te));

for i=1:length(Tn)
    for j=1:length(Te)
        
        Ln=Tn(i)/100.0;
        Le=Te(j)/100.0;
        
        Lc=1-Le-Ln;
        
        L = 1e6*(Le*exp(the*tt)+Lc*exp(thc*tt)+Ln*exp(thn*tt));
        
        inp1=find(L<1e6/70);
        inpc=find(L<400);
        inH1=find(L<500);
        inHc=find(L<100);

        p1(i,j)=t(inp1(1));
        pc(i,j)=t(inpc(1));
        H1(i,j)=t(inH1(1));
        Hc(i,j)=t(inHc(1));
        
    end
end

% zoomed calculation

Tnz=linspace(0,10,100); %percent naive
Tez=linspace(0,40,100); %percent em
 
t=linspace(0,10,10000);

tt=t*365; %time in days

for i=1:length(Tnz)
    for j=1:length(Tez)
        
        Ln=Tnz(i)/100.0;
        Le=Tez(j)/100.0;
        
        Lc=1-Le-Ln;
        
        L = 1e6*(Le*exp(the*tt)+Lc*exp(thc*tt)+Ln*exp(thn*tt));
        
        inp1=find(L<1e6/70);
        inpc=find(L<400);
        inH1=find(L<500);
        inHc=find(L<100);

        p1z(i,j)=t(inp1(1));
        pcz(i,j)=t(inpc(1));
        H1z(i,j)=t(inH1(1));
        Hcz(i,j)=t(inHc(1));
        
    end
end
%%
figuresize(7,4,'inches')

subplot(231)
[~,hC]=contourf(Te,Tn,p1,25);
shading interp
set(hC,'LineStyle','none');
colorbar
caxis([1,7])

subplot(232)
[~,hC]=contourf(Te,Tn,H1,25);
shading interp
set(hC,'LineStyle','none');
colorbar
caxis([1,7])

% subplot(243)
% [~,hC]=contourf(Te,Tn,pc,25);
% shading interp
% set(hC,'LineStyle','none');
% colorbar
% caxis([1,7])

subplot(233)
[~,hC]=contourf(Te,Tn,Hc,25);
shading interp
set(hC,'LineStyle','none');
colorbar
caxis([1,7])

subplot(234)
[~,hC]=contourf(Tez,Tnz,p1z,10);
shading interp
%set(hC,'LineStyle','none');
colorbar
caxis([0,2])

subplot(235)
[~,hC]=contourf(Tez,Tnz,H1z,10);
shading interp
%set(hC,'LineStyle','none');
xlim([0,40]); 
set(gca,'XTick',linspace(0,40,5))
ylim([0,40])
set(gca,'YTick',linspace(0,40,5))
colorbar
caxis([2,5])

% subplot(247)
% [~,hC]=contourf(Tez,Tnz,pcz);
% shading interp
% %set(hC,'LineStyle','none');
% colorbar
% caxis([2,5])

subplot(236)
[~,hC]=contourf(Tez,Tnz,Hcz,10);
shading interp
%set(hC,'LineStyle','none');
xlim([0,40]); 
set(gca,'XTick',linspace(0,40,5))
ylim([0,40])
set(gca,'YTick',linspace(0,40,5))
colorbar
caxis([2,5])

colormap('jet')

cs = findall(gcf); for k = 1:numel(cs); try set(cs(k), 'FontSize',14); end; end


print('PIC_subsets','-dtiff','-r600')


%% plot for florian

%% figure 3 plots of heterogeneous latent pools

thL=-5.2e-4;
xi = 5.7e-5;

aTn = 0.002;
dTn = aTn-thL-xi;

aTc = 0.015;
dTc = aTc-thL-xi;

aTe = 0.047;
dTe = aTe-thL-xi;

thc = aTc/10-dTc-xi;
the = aTe/10-dTe-xi;
thn = aTn/10-dTn-xi;

%full calculation

Tn=linspace(0,40,100); %percent naive
Te=linspace(0,40,100); %percent em
 
t=linspace(0,10,10000);

tt=t*365; %time in days

L=zeros(length(Tn),length(Te));

for i=1:length(Tn)
    for j=1:length(Te)
        
        Ln=Tn(i)/100.0;
        Le=Te(j)/100.0;
        
        Lc=1-Le-Ln;
        
        L = 1e6*(Le*exp(the*tt)+Lc*exp(thc*tt)+Ln*exp(thn*tt));
        
        inp1=find(L<1e6/70);
        inpc=find(L<400);
        inH1=find(L<500);
        inHc=find(L<100);

        p1(i,j)=t(inp1(1));
        pc(i,j)=t(inpc(1));
        H1(i,j)=t(inH1(1));
        Hc(i,j)=t(inHc(1));
        
    end
end

%% Fig 3 actually plot contour plots

figure()
figuresize(5,1.5,'inches')

subplot(131)
[~,hC]=contourf(Te,Tn,p1,25);
shading interp
set(hC,'LineStyle','none');
xlim([0,40]); 
set(gca,'XTick',linspace(0,40,5))
ylim([0,40])
set(gca,'YTick',linspace(0,40,5))
%colorbar
caxis([1,7])

subplot(132)
[~,hC]=contourf(Te,Tn,H1,25);
shading interp
set(hC,'LineStyle','none');
xlim([0,40]); 
set(gca,'XTick',linspace(0,40,5))
ylim([0,40])
set(gca,'YTick',linspace(0,40,5))
%colorbar
caxis([1,7])

subplot(133)
[~,hC]=contourf(Te,Tn,Hc,25);
shading interp
set(hC,'LineStyle','none');
xlim([0,40]); 
set(gca,'XTick',linspace(0,40,5))
ylim([0,40])
set(gca,'YTick',linspace(0,40,5))
%colorbar
caxis([1,7])

colormap('jet')
cs = findall(gcf); for k = 1:numel(cs); try set(cs(k), 'FontSize',10); end; end

print('PIC_subsets','-dtiff','-r600')

%%
figure()
figuresize(4,5,'inches')

LFe=0.4; LFe_t=1e6*(LFe*exp(the*tt));
LFc=0.59; LFc_t=1e6*(LFc*exp(thc*tt));
LFn=1-LFe-LFc; LFn_t=1e6*(LFn*exp(thn*tt));


subplot(211)
semilogy(t,LFe_t,'color',[31, 119, 180]/255,'linewidth',2)
hold on
semilogy(t,LFc_t,'color',[152, 223, 138]/255,'linewidth',2)
semilogy(t,LFn_t,'color',[214, 39, 40]/255,'linewidth',2)
semilogy(t,LFe_t+LFc_t+LFn_t,'color','k','linewidth',3,'linestyle',':')
semilogy(t,ones(length(t),1)*1.4e4,'color','black','linestyle','-')
semilogy(t,ones(length(t),1)*500,'color','black','linestyle','-')
semilogy(t,ones(length(t),1)*100,'color','black','linestyle','-')
hold off
legend('T_{em} (40%)','T_{cm} (59%)','T_{n }  (1%)',' sum')
xlim([0,6])
ylim([80,2e6])
set(gca,'YTick',logspace(2,8,7))
set(gca,'XTick',0:6)

LFe=0.4; LFe_t=1e6*(LFe*exp(the*tt));
LFc=0.5; LFc_t=1e6*(LFc*exp(thc*tt));
LFn=1-LFe-LFc; LFn_t=1e6*(LFn*exp(thn*tt));

subplot(212)
semilogy(t,LFe_t,'color',[31, 119, 180]/255,'linewidth',2)
hold on
semilogy(t,LFc_t,'color',[152, 223, 138]/255,'linewidth',2)
semilogy(t,LFn_t,'color',[214, 39, 40]/255,'linewidth',2)
semilogy(t,LFe_t+LFc_t+LFn_t,'color','k','linewidth',3,'linestyle',':')
semilogy(t,ones(length(t),1)*1.4e4,'color','black','linestyle','-')
semilogy(t,ones(length(t),1)*500,'color','black','linestyle','-')
semilogy(t,ones(length(t),1)*100,'color','black','linestyle','-')
hold off
xlim([0,6])
ylim([80,1e6])
legend('T_{em} (40%)','T_{cm} (50%)','T_{n }  (10%)',' sum')
set(gca,'XTick',0:6)
set(gca,'YTick',logspace(2,8,7))

cs = findall(gcf); for k = 1:numel(cs); try set(cs(k), 'FontSize',14); end; end

print('PIC_subsets_lines','-dpdf','-r600')

%% simulate the waning ap efficacy
clf

tspan=linspace(0,10*365,1e4); %10 years
    
fc=0.55;
fn=0.01;
fe=1-fc-fn;

L0=1e6;

thL=-5.2e-4;
xi = 1.1e-5;

aTn = 0.002;
dTn = aTn-thL-xi;

aTc = 0.015;
dTc = aTc-thL-xi;

aTe = 0.047;
dTe = aTe-thL-xi;

ep0=5;
ap_decay = [0 0.01/30 0.05/30 0.1/30 0.2/30];

figuresize(5,3,'inches')
hold on
for i =1:length(ap_decay)
    
    phi=ap_decay(i);
    [t,ye] = ode45(@(t,y) (aTe/(1+(ep0-1)*exp(-phi*t))-dTe-xi)*y, tspan, L0);
    [t,yc] = ode45(@(t,y) (aTc/(1+(ep0-1)*exp(-phi*t))-dTc-xi)*y, tspan, L0);
    [t,yn] = ode45(@(t,y) (aTn/(1+(ep0-1)*exp(-phi*t))-dTn-xi)*y, tspan, L0);

    semilogy(t/365, ye*fe + yc*fc + yn*fn,'Linewidth',3)

end

line(get(gca,'XLim'),[100 100],'Color','k')
line(get(gca,'XLim'),[500 500],'Color','k')
line(get(gca,'XLim'),[1.4e4 1.4e4],'Color','k')
hold off
set(gca,'yScale','log')
ylim([0,1e6])
set(gca,'xTick',0:10)
set(gca,'yTick',logspace(1,6,6))
legend('no decrease','1% per mo.','5% per mo.','10% per mo.','20% per mo.')

print('potency_decrease','-dpdf','-r600')
