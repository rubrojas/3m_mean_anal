
%_________________________________________________________________________%

% SINGLE DAY PLOTTER
% ALL NOTES FROM COMPILE_GAIN_PLOTS
%compile and plot gains all days
%It plots Bs,bphi2, for the days indicated. It alsos pauses so you
%can go checking them day by day.
%06/22/2020 MODIFIED TO INCLUDE ONLY m0 since according to plots in
%plot_garray these are the main ones in amplitude.
% function single_day_amp_plot(run)
%_________________________________________________________________________%
%CHOOSE YOUR DAY WORKSPACE IN HERE
% run = '081722_1';
day = [run 'amp.mat']; %when using function

%_________________________________________________________________________%
%CHOOSE NUMBER OF FIGURES

fignum1 = 14;
fignum2 = 15;
    
    
%_________________________________________________________________________%
%CHOOSE IF WANNA PLOT ALL GCOEFF
GAUSS = 1;

%_________________________________________________________________________%
%CHOOSE TO SAVE PLOTS IN FOLDER 
SALVA = 0;

%_________________________________________________________________________%
%size of PLOTS
pz = 1.5;

%FONT SIZE
fz = 14;

%scatter Size

sz = 50;
%_________________________________________________________________________%
%LOADING 
cd(['/Users/rubenrojas/Documents/3M/Amplification/3m_folder_files/' day(1:6)]);
masterdir=dir([pwd '/' day]);
load(day, '-regexp', '^(?!record|masterdir)\w') %LOADING EVERYTHING EXCEPT THE FKING HEAVY RECORD MATRIX
%_________________________________________________________________________%
%SELECT THE GAUSS COEFF TO PLOT IN THE FINGER PLOT
bla=bl1m0;
blb=bl1m1;
blc=bl3m0;
bld=bl2m1;
ble=bl4m0;
blf=bl3m2;
blg=bl4m2;

%Var names for legend
% na=getname(b
bla_var='$B_{1}^0$';
blb_var='$B_{1}^1$';
blc_var='$B_{3}^0$';
bld_var='$B_{2}^1$';
ble_var='$B_{4}^0$';
blf_var='$B_{3}^2$';
blg_var='$B_{4}^2$';

%_________________________________________________________________________%
%THIS IS WHEN ALPHA DATA ARRIVES
% day=str2double(folder);
% if day(end-2:end) < 18
% masterdir=dir('/Users/rubenrojas/Documents/3M/Amplification/3m_smooth/records');
% else
% masterdir=dir('/Users/rubenrojas/Documents/3M/Amplification/3m_alpha/');
% end


% set(gcf,'position',[0 0 600 800]);%left bottom width height
cmap=[1,0,0;1,0.396000000000000,0;1,0.792000000000000,0;0.812000000000000,1,0;0.416000000000000,1,0;0.0200000000000000,1,0;0,1,0.376000000000000;0,1,0.772000000000000;0,0.832000000000000,1;0,0.436000000000000,1;0,0.0400000000000000,1;0.356000000000000,0,1;0.752000000000001,0,1;1,0,0.852000000000000;1,0,0.456000000000000;1,0,0.0600000000000005];
% cmap9=[1,0,0;1,0.400000000000000,0;1,0.800000000000000,0;0.500000000000000,1,0;0,1,0.772000000000000;0,0.832000000000000,1;0,0.0400000000000000,1;0.752000000000001,0,1;1,0,0.456000000000000];
% cmap5=[1,0,0;1,0.800000000000000,0;0.500000000000000,1,0;0,0.0400000000000000,1;0.752000000000001,0,1];
cmaplab=[0, 0.4470, 0.7410;0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250;0.4940, 0.1840, 0.556;0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.933;0.6350, 0.0780, 0.1840];



% % cleaning spoureous data ro1=-1 or ro1=0
%         a=ro1==0;b=ro1==-1;c=bext<=10;
%         ro1(a)=NaN; ro1(b)=NaN;ro1(c)=NaN;
%         bext(a)=NaN; bext(b)=NaN;bext(c)=NaN;
%         rm(a)=NaN; rm(b)=NaN;rm(c)=NaN;

%identifying id Cuad, Dipole, or Single coil
if str2double(folder(5:6)) >= 15  && str2double(folder(5:6)) < 21
    var_label='D';
    if strcmp(folder,'032416') || strcmp(folder,'111715')
        var_label='C';
    else
    end
else
    var_label='S';
end

if str2double(folder(5:6)) >= 21
    var_label='D';
end




%CREATING GARRAY----------------------------------------------------------
%
%         %             m=0      m=1(cos),m=-1(sin)      m=2(cos),m=-2(sin)      m=3(cos),m=-3(sin)      m=4(cos),m=-4(sin)
%         gcarray= {bla,{bl1mc1,bl1ms1},       {NaN,NaN},              {NaN,NaN},              {NaN,NaN};        %l=1
%             blb,{bl2mc1,bl2ms1},{bl2mc2,bl2ms2},       {NaN,NaN},              {NaN,NaN};        %l=2
%             blc,{bl3mc1,bl3ms1},{bl3mc2,bl3ms2},{bl3mc3,bl3ms3},       {NaN,NaN};        %l=3
%             bld,{bl4mc1,bl4ms1},{bl4mc2,bl4ms2},{bl4mc3,bl4ms3},{bl4mc4,bl4ms4}};%l=4
%
%THIS IS ANOTHER GARRAY USING DAN NORM cos2 + sin2 for the glm
%          m=0   m=1   m=2   m=3   m=4
gcarray= {bl1m0,bl1m1,NaN  ,NaN  ,NaN;        %l=1
          bl2m0,bl2m1,bl2m2,NaN  ,NaN;        %l=2
          bl3m0,bl3m1,bl3m2,bl3m3,NaN;        %l=3
          bl4m0,bl4m1,bl4m2,bl4m3,bl4m4};     %l=4


%% PLOTTING ALL GCAOEFF TO COMPARE----------------------------------------------------------------
if GAUSS ==1
figure(fignum2);clf;
set(gcf,'color','w');


for i=1:4
    
    subplot(3,4,i);
    
    for j=1:i+1
        
        plot(ro1,gcarray{i,j},'-','color',cmaplab(j,:),'linewidth',pz);hold on;
        scatter(ro1,gcarray{i,j},100,'.','markeredgecolor',cmaplab(j,:),'handlevisibility','off');
    end
    
    xlabel('$Ro^{-1}$','FontName','Times','Interpreter','latex');
    if norma == 1; ylabel('$B_{l}^m/B_0$','FontName','Times','Interpreter','latex');else;ylabel('$B_{l}^m (\mathrm{G})$','FontName','Times','Interpreter','latex');end
    title([ day(1:end-9) '.' day(8) '(' var_label ')' ' l = ' num2str(i)],'Interpreter','latex');
    if i==4
        n=legend({'$m=0$','$m=1$','$m=2$','$m=3$','$m=4$'},'Interpreter','latex');
        set(n,'position',get(n, 'position').*[1.11 1 1 1]);
    end
    box on
    hold off
    
    subplot(3,4,i+4);
    
    for j=1:i+1
        
        plot(rm,gcarray{i,j},'-','color',cmaplab(j,:),'linewidth',pz);hold on;
        scatter(rm,gcarray{i,j},100,'.','markeredgecolor',cmaplab(j,:),'handlevisibility','off');
    end
    
    
    if norma == 1; ylabel('$B_{l}^m/B_0$','FontName','Times','Interpreter','latex');else;ylabel('$B_{l}^m (\mathrm{G})$','FontName','Times','Interpreter','latex');end
    xlabel('$Rm$','FontName','Times','Interpreter','latex');
    box on
    hold off
    
    subplot(3,4,i+8);
    
    for j=1:i+1
        
        plot(bext,gcarray{i,j},'-','color',cmaplab(j,:),'linewidth',pz);hold on;
        scatter(bext,gcarray{i,j},100,'.','markeredgecolor',cmaplab(j,:),'handlevisibility','off');
        
    end
    if norma == 1; ylabel('$B_{l}^m/B_0$','FontName','Times','Interpreter','latex');else;ylabel('$B_{l}^m (\mathrm{G})$','FontName','Times','Interpreter','latex');end
    xlabel('$B_{ext}$','FontName','Times','Interpreter','latex');
    
    box on
    hold off
    
    
end
%     if SALVA == 1
%     saveas(gcf,['/Users/rubenrojas/Documents/3M/Amplification/3m_folder_files/' masterdir.name(1:end-9) '/' masterdir.name(1:end-7) 'gauss.png'],'png');
%     end

else
end

%% THIS IS FOR PLOT FINGER --------------------------------
figure(fignum1);clf;
set(gcf,'color','w');

subplot(3,3,1);hold on;

v=scatter(ro1,bphi,sz,'o','Markeredgecolor','r'); hold on;
l=scatter(ro1,brad,sz,'*','Markeredgecolor','r'); hold on;
p=plot(ro1,brad,'color','r','linewidth',pz);
plot(ro1,bphi,'color','r','linewidth',pz);
h=scatter(ro1(find(~isnan(ro1),1,'first')),brad(find(~isnan(ro1),1,'first')),sz,'O','filled','Markerfacecolor','k');hold on;
scatter(ro1(find(~isnan(ro1),1,'first')),bphi(find(~isnan(ro1),1,'first')),sz,'O','filled','Markerfacecolor','k');

if str2double(folder(5:6)) >= 21
scatter(ro1,bphi2,sz,'o','Markeredgecolor','k');
scatter(ro1,brad2,sz,'*','Markeredgecolor','k');
j=plot(ro1,brad2,'color','k','linewidth',pz);
plot(ro1,bphi2,'color','k','linewidth',pz);
scatter(ro1(find(~isnan(ro1),1,'first')),brad2(find(~isnan(ro1),1,'first')),sz,'O','filled','Markerfacecolor','k');hold on;
scatter(ro1(find(~isnan(ro1),1,'first')),bphi2(find(~isnan(ro1),1,'first')),sz,'O','filled','Markerfacecolor','k');
end

%axis([-0.06 -0.05 -1 2])
%hold off
xlabel('$Ro^{-1}$','FontName','Times','Interpreter','latex');
if norma == 1;ylabel('$B_{i}/B_0$','FontName','Times','Interpreter','latex');else;ylabel('$B_i (\mathrm{G})$','FontName','Times','Interpreter','latex');end
title([day(1:end-9) '.' day(8) '(' var_label ')'],'Interpreter','latex');
if str2double(folder(5:6)) >= 21
m=legend([l,v,p,j,h],{'$B_r$','$B_{\varphi}$','Port D','Port B','$t_0$'},'location','northeast','Interpreter','latex');
else
m=legend([l,v,h],{'$B_r$','$B_{\varphi}$','$t_0$'},'location','northeast','Interpreter','latex');
   
end
set(m,'position',[0.011    0.7812    0.0655    0.0935]);
set(gca,'linewidth',pz,'FontName','Times New Roman','FontSize',fz)
box on
hold off

subplot(3,3,4);hold on;
scatter(rm,brad,sz,'*','Markeredgecolor','r');hold on;
scatter(rm,bphi,sz,'o','Markeredgecolor','r');
plot(rm,brad,'color','r','linewidth',pz);
plot(rm,bphi,'color','r','linewidth',pz);
scatter(rm(find(~isnan(rm),1,'first')),brad(find(~isnan(rm),1,'first')),sz,'O','filled','Markerfacecolor','k');hold on;
scatter(rm(find(~isnan(rm),1,'first')),bphi(find(~isnan(rm),1,'first')),sz,'O','filled','Markerfacecolor','k');

if str2double(folder(5:6)) >= 21
scatter(rm,brad2,sz,'*','Markeredgecolor','k');hold on;
scatter(rm,bphi2,sz,'o','Markeredgecolor','k');
plot(rm,brad2,'color','k','linewidth',pz);
plot(rm,bphi2,'color','k','linewidth',pz);
scatter(rm(find(~isnan(rm),1,'first')),brad2(find(~isnan(rm),1,'first')),sz,'O','filled','Markerfacecolor','k');hold on;
scatter(rm(find(~isnan(rm),1,'first')),bphi2(find(~isnan(rm),1,'first')),sz,'O','filled','Markerfacecolor','k');
end

if norma == 1;ylabel('$B_{i}/B_0$','FontName','Times','Interpreter','latex');else;ylabel('$B_i (\mathrm{G})$','FontName','Times','Interpreter','latex');end
xlabel('$Rm$','FontName','Times','Interpreter','latex');
set(gca,'linewidth',pz,'FontName','Times New Roman','FontSize',fz)
title([day(1:end-9) '.' day(8) '(' var_label ')'],'Interpreter','latex');
box on
hold off

subplot(3,3,7);hold on;
scatter(bext,brad,sz,'*','Markeredgecolor','r'); hold on;
scatter(bext,bphi,sz,'o','Markeredgecolor','r');
plot(bext,brad,'color','r','linewidth',pz);
plot(bext,bphi,'color','r','linewidth',pz);
scatter(bext(find(~isnan(bext),1,'first')),brad(find(~isnan(bext),1,'first')),sz,'O','filled','Markerfacecolor','k');hold on;
scatter(bext(find(~isnan(bext),1,'first')),bphi(find(~isnan(bext),1,'first')),sz,'O','filled','Markerfacecolor','k');

if str2double(folder(5:6)) >= 21
scatter(bext,brad2,sz,'*','Markeredgecolor','k'); hold on;
scatter(bext,bphi2,sz,'o','Markeredgecolor','k');
plot(bext,brad2,'color','k','linewidth',pz);
plot(bext,bphi2,'color','k','linewidth',pz);
scatter(bext(find(~isnan(bext),1,'first')),brad2(find(~isnan(bext),1,'first')),sz,'O','filled','Markerfacecolor','k');hold on;
scatter(bext(find(~isnan(bext),1,'first')),bphi2(find(~isnan(bext),1,'first')),sz,'O','filled','Markerfacecolor','k');
end

%scatter(tim,ro1,sz,'*','Markeredgecolor','k');
%scatter(tim,ro1,sz,'o','Markeredgecolor','k');
%xlim([0 max(bext)])
if norma == 1;ylabel('$B_{i}/B_0$','FontName','Times','Interpreter','latex');else;ylabel('$B_i (\mathrm{G})$','FontName','Times','Interpreter','latex');end
xlabel('$B_{ext}$','FontName','Times','Interpreter','latex');
set(gca,'linewidth',pz,'FontName','Times New Roman','FontSize',fz)
title([day(1:end-9) '.' day(8) '(' var_label ')'],'Interpreter','latex');

box on
hold off
%     end

%THIS IS FOR TORQUE -----------------------------------------------------
subplot(3,3,2);hold on;
%scatter3(ro1,bphi2,ones(1,length(ro1))*str2num(folder),sz,'o','Markeredgecolor','k'); %here I normalized
scatter(ro1,gginf,sz,'O','filled','Markerfacecolor','b'); hold on;
%         plot(ro1,gginf,'color','k');
%axis([-0.06 -0.05 -1 2])
%         hold off
plot(ro1,gginf,'color','b','linewidth',pz);
scatter(ro1(find(~isnan(ro1),1,'first')),gginf(find(~isnan(ro1),1,'first')),sz,'O','filled','Markerfacecolor','k');

xlabel('$Ro^{-1}$','FontName','Times','Interpreter','latex');
ylabel('$G/G_{\infty}$','FontName','Times','Interpreter','latex');
% ylabel('$G$','FontName','Times','Interpreter','latex');

title([day(1:end-9) '.' day(8) '(' var_label ')'],'Interpreter','latex');
set(gca,'linewidth',pz,'FontName','Times New Roman','FontSize',fz)
box on
hold off

subplot(3,3,5);hold on;
%scatter3(rm,brad2,ones(1,length(ro1))*str2num(folder),sz,'*','Markeredgecolor','k');%here I normalized
%scatter3(rm,bphi2,ones(1,length(ro1))*str2num(folder),sz,'o','Markeredgecolor','k'); %here I normalized
scatter(rm,gginf,sz,'O','filled','Markerfacecolor','b'); hold on;
plot(rm,gginf,'color','b','linewidth',pz);
scatter(rm(find(~isnan(rm),1,'first')),gginf(find(~isnan(rm),1,'first')),sz,'O','filled','Markerfacecolor','k');

ylabel('$G/G_{\infty}$','FontName','Times','Interpreter','latex');
% ylabel('$G$','FontName','Times','Interpreter','latex');

xlabel('$Rm$','FontName','Times','Interpreter','latex');
set(gca,'linewidth',pz,'FontName','Times New Roman','FontSize',fz)
box on
hold off

subplot(3,3,8);hold on;
scatter(bext,gginf,sz,'O','filled','Markerfacecolor','b'); hold on;
ylabel('$G/G_{\infty}$','FontName','Times','Interpreter','latex');
% ylabel('$G$','FontName','Times','Interpreter','latex');
xlabel('$B_{ext}$','FontName','Times','Interpreter','latex');
plot(bext,gginf,'color','b','linewidth',pz);
scatter(bext(find(~isnan(bext),1,'first')),gginf(find(~isnan(bext),1,'first')),sz,'O','filled','Markerfacecolor','k');
set(gca,'linewidth',pz,'FontName','Times New Roman','FontSize',fz)
title([day(1:end-9) '.' day(8) '(' var_label ')'],'Interpreter','latex');
box on    
hold off



%THIS IS FOR SELECTED GAUSS COEFF ------------------------------------------
subplot(3,3,3);hold on;

%scatter3(ro1,bphi2,ones(1,length(ro1))*str2num(folder),sz,'o','Markeredgecolor','k'); %here I normalized
scatter(ro1,bla,sz,'filled','Markerfacecolor',cmap(1,:)); hold on;
scatter(ro1,blb,sz,'filled','Markerfacecolor',cmap(3,:));
scatter(ro1,blc,sz,'filled','Markerfacecolor',cmap(5,:));
scatter(ro1,bld,sz,'filled','Markerfacecolor',cmap(11,:));
scatter(ro1,ble,sz,'filled','Markerfacecolor',cmap(2,:));
scatter(ro1,blf,sz,'filled','Markerfacecolor',cmap(8,:));
scatter(ro1,blg,sz,'filled','Markerfacecolor',cmap(15,:));
a=plot(ro1,bla,'color',cmap(1,:),'linewidth',pz);hold on;
b=plot(ro1,blb,'color',cmap(3,:),'linewidth',pz);
c=plot(ro1,blc,'color',cmap(5,:),'linewidth',pz);
d=plot(ro1,bld,'color',cmap(11,:),'linewidth',pz);
e=plot(ro1,ble,'color',cmap(2,:),'linewidth',pz);
f=plot(ro1,blf,'color',cmap(8,:),'linewidth',pz);
h=plot(ro1,blg,'color',cmap(15,:),'linewidth',pz);
scatter(ro1(find(~isnan(ro1),1,'first')),bla(find(~isnan(ro1),1,'first')),'Markerfacecolor','k');hold on;
scatter(ro1(find(~isnan(ro1),1,'first')),blb(find(~isnan(ro1),1,'first')),'Markerfacecolor','k');
scatter(ro1(find(~isnan(ro1),1,'first')),blc(find(~isnan(ro1),1,'first')),'Markerfacecolor','k');
scatter(ro1(find(~isnan(ro1),1,'first')),bld(find(~isnan(ro1),1,'first')),'Markerfacecolor','k');
scatter(ro1(find(~isnan(ro1),1,'first')),ble(find(~isnan(ro1),1,'first')),'Markerfacecolor','k');
scatter(ro1(find(~isnan(ro1),1,'first')),blf(find(~isnan(ro1),1,'first')),'Markerfacecolor','k');
scatter(ro1(find(~isnan(ro1),1,'first')),blg(find(~isnan(ro1),1,'first')),'Markerfacecolor','k');
xlabel('$Ro^{-1}$','FontName','Times','Interpreter','latex');
if norma == 1; ylabel('$B_{l}^m/B_0$','FontName','Times','Interpreter','latex');else;ylabel('$B_{l}^m (\mathrm{G})$','FontName','Times','Interpreter','latex');end
    
title([day(1:end-9) '.' day(8) '(' var_label ')'],'Interpreter','latex');
n=legend([a,b,c,d,e,f,h],{bla_var,blb_var,blc_var,bld_var,ble_var,blf_var,blg_var},'Interpreter','latex');
% n=legend([a,b,c,d,e,f,g],{'$B_{1}^0$','$B_{1}^1$','$B_{3}^0$','$B_{4}^0$', '$B_{2}^1$','$B_{3}^2$','$B_{4}^3$'},'Interpreter','latex');
set(n,'position',get(n, 'position').*[1.1 1 1 1]);
set(gca,'linewidth',pz,'FontName','Times New Roman','FontSize',fz)
box on

hold off

subplot(3,3,6);hold on;
scatter(rm,bla,sz,'filled','Markerfacecolor',cmap(1,:)); hold on;
scatter(rm,blb,sz,'filled','Markerfacecolor',cmap(3,:));
scatter(rm,blc,sz,'filled','Markerfacecolor',cmap(5,:));
scatter(rm,bld,sz,'filled','Markerfacecolor',cmap(11,:));
scatter(rm,ble,sz,'filled','Markerfacecolor',cmap(2,:));
scatter(rm,blf,sz,'filled','Markerfacecolor',cmap(8,:));
scatter(rm,blg,sz,'filled','Markerfacecolor',cmap(15,:));

plot(rm,bla,'color',cmap(1,:),'linewidth',pz);hold on;
plot(rm,blb,'color',cmap(3,:),'linewidth',pz);
plot(rm,blc,'color',cmap(5,:),'linewidth',pz);
plot(rm,bld,'color',cmap(11,:),'linewidth',pz);
plot(rm,ble,'color',cmap(2,:),'linewidth',pz);
plot(rm,blf,'color',cmap(8,:),'linewidth',pz);
plot(rm,blg,'color',cmap(15,:),'linewidth',pz);
scatter(rm(find(~isnan(rm),1,'first')),bla(find(~isnan(rm),1,'first')),'Markerfacecolor','k');hold on;
scatter(rm(find(~isnan(rm),1,'first')),blb(find(~isnan(rm),1,'first')),'Markerfacecolor','k');
scatter(rm(find(~isnan(rm),1,'first')),blc(find(~isnan(rm),1,'first')),'Markerfacecolor','k');
scatter(rm(find(~isnan(rm),1,'first')),bld(find(~isnan(rm),1,'first')),'Markerfacecolor','k');
scatter(rm(find(~isnan(rm),1,'first')),ble(find(~isnan(rm),1,'first')),'Markerfacecolor','k');
scatter(rm(find(~isnan(rm),1,'first')),blf(find(~isnan(rm),1,'first')),'Markerfacecolor','k');
scatter(rm(find(~isnan(rm),1,'first')),blg(find(~isnan(rm),1,'first')),'Markerfacecolor','k');
if norma == 1; ylabel('$B_{l}^m/B_0$','FontName','Times','Interpreter','latex');else;ylabel('$B_{l}^m (\mathrm{G})$','FontName','Times','Interpreter','latex');end
xlabel('$Rm$','FontName','Times','Interpreter','latex');
set(gca,'linewidth',pz,'FontName','Times New Roman','FontSize',fz)
title([day(1:end-9) '.' day(8) '(' var_label ')'],'Interpreter','latex');

box on
hold off

subplot(3,3,9);
scatter(bext,bla,sz,'filled','Markerfacecolor',cmap(1,:)); hold on;
scatter(bext,blb,sz,'filled','Markerfacecolor',cmap(3,:));
scatter(bext,blc,sz,'filled','Markerfacecolor',cmap(5,:));
scatter(bext,bld,sz,'filled','Markerfacecolor',cmap(11,:));
scatter(bext,ble,sz,'filled','Markerfacecolor',cmap(2,:));
scatter(bext,blf,sz,'filled','Markerfacecolor',cmap(8,:));
scatter(bext,blg,sz,'filled','Markerfacecolor',cmap(15,:));

plot(bext,bla,'color',cmap(1,:),'linewidth',pz);hold on;
plot(bext,blb,'color',cmap(3,:),'linewidth',pz);
plot(bext,blc,'color',cmap(5,:),'linewidth',pz);
plot(bext,bld,'color',cmap(11,:),'linewidth',pz);
plot(bext,ble,'color',cmap(2,:),'linewidth',pz);
plot(bext,blf,'color',cmap(8,:),'linewidth',pz);
plot(bext,blg,'color',cmap(15,:),'linewidth',pz);

scatter(bext(find(~isnan(bext),1,'first')),bla(find(~isnan(bext),1,'first')),'Markerfacecolor','k');hold on;
scatter(bext(find(~isnan(bext),1,'first')),blb(find(~isnan(bext),1,'first')),'Markerfacecolor','k');
scatter(bext(find(~isnan(bext),1,'first')),blc(find(~isnan(bext),1,'first')),'Markerfacecolor','k');
scatter(bext(find(~isnan(bext),1,'first')),bld(find(~isnan(bext),1,'first')),'Markerfacecolor','k');
scatter(bext(find(~isnan(bext),1,'first')),ble(find(~isnan(bext),1,'first')),'Markerfacecolor','k');
scatter(bext(find(~isnan(bext),1,'first')),blf(find(~isnan(bext),1,'first')),'Markerfacecolor','k');
scatter(bext(find(~isnan(bext),1,'first')),blg(find(~isnan(bext),1,'first')),'Markerfacecolor','k');


if norma == 1; ylabel('$B_{l}^m/B_0$','FontName','Times','Interpreter','latex');else;ylabel('$B_{l}^m (\mathrm{G})$','FontName','Times','Interpreter','latex');end
xlabel('$B_{ext}$','FontName','Times','Interpreter','latex');
set(gca,'linewidth',pz,'FontName','Times New Roman','FontSize',fz)
title([day(1:end-9) '.' day(8) '(' var_label ')'],'Interpreter','latex');

box on
hold off




%THIS IS FOR SEEING AND SAVING PLOTS
if SALVA == 1

%        saveas(gcf,['/Users/rubenrojas/Documents/3M/Amplification/3m_smooth/figures_amp/' masterdir.name(1:end-7) 'finger.png'],'png');
%    pause;
       saveas(gcf,['/Users/rubenrojas/Documents/3M/Amplification/3m_rough/figures_amp/' masterdir.name(1:end-7) 'finger.png'],'png');

%        saveas(gcf,['/Users/rubenrojas/Documents/3M/Amplification/3m_folder_files/' masterdir.name(1:end-9) '/' masterdir.name(1:end-7) 'finger.png'],'png');
end

clearvars a b c d e f n

%     end


