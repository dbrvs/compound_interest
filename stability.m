% stability analysis of the SLAV differential equations
% DBR 12/15

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

ep   = 0.9;   %effectiveness of ART from 0->1 1 being 100% effective []
betaE=beta*(1-ep);

JeqVFE =[[-dS 0 0 -betaE*aS/dS];...
         [0 thL 0 tau*betaE*aS/dS];...
         [0 xi aA-dA (1-tau)*betaE*aS/dS];...
         [0 0 p -g]];
     
for i=1:99
    
ep(i)=i/100;
betaE=beta*(1-ep(i));
    
S=g*dA/betaE/p/fL;
L=tau/thL*(g*dA*dS/betaE/p/fL-aS);
A=aS*fL/dA-g*dS/p/betaE;
V=p*aS*fL/g/dA-dS/betaE;

Jeq =[[-dS-betaE*V 0 0 -betaE*S];...
      [tau*betaE*V thL 0 tau*betaE*S];...
      [(1-tau)*betaE*V xi aA-dA (1-tau)*betaE*S];...
      [0 0 p -g]];

eD(:,i)=eig(Jeq);

slav(:,i)=[S L A V];

%ctime=1/min(abs(eD))/365;

end

%%
%    legend('Seq','Leq','Aeq','Veq','Location','SouthWest')
names = ['S*'; 'L*'; 'A*' ;'V*'];

ylimz=[[-30,5];[-1,1];[-1,0.5];[-6e4,100]];

figure(1)
figuresize(6,7,'inches')
for j=1:4
    subplot(4,2,j*2-1)
    plot(ep,slav(j,:),'linewidth',3)
    hold on
    line([cep cep],get(gca,'YLim'),'Color','k','LineStyle','--')
    hold off
    ylim([0,max(slav(j,:))])
    set(gca,'xTick',linspace(0,1,5))
    grid on
    ylabel(names(j,:), 'Fontsize',12)
    if j==4
        xlabel('ART efficacy \epsilon^{ART}')
    end
 
    subplot(4,2,j*2)
    plot(ep,real(eD(j,:)),'linewidth',3)
    hold on
    plot(ep,imag(eD(j,:)),'linewidth',3)
    line([cep cep],get(gca,'YLim'),'Color','k','LineStyle','--')
    set(gca,'xTick',linspace(0,1,5))
    hold off
    if j==4
        legend('real','imag','Location','SouthWest')
        xlabel('ART efficacy \epsilon^{ART}')
    end
    ylabel(['\lambda_' +num2str(j) ])
    %ylim(ylimz(j,:))
    grid on
    
end

