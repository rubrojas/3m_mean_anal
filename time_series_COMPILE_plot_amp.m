% % % time series plot for amplification compiling
%Changing this to do many plots on top for a particular run. APRIL21 Ruben
%PLOTTING RESULTS FOR MANY POINTS ON A SINGLE RUN, TO COMPARE SHORT TERM DYNAMICS AMONG THEM

%Preeliminares for run day------------------------------------------------------
if exist('run','var'); passrun=run;else;  end %prelocating the old 'run' to make sure you loaded the record matrix

% run = '020922_2'; %run day

%this is for selecting the point within each run
% path = [2,4,5,6,7,8];
path = 'full';

%type of ramp: CHOOSE TYPE OF RAMP ------------------------------------
% ramp = 'ros'; %roosby ramp
ramp ='rey'; %reynolds ramp
% ramp = 'mag'; %mag ramp

%chooses the porcentaje of the time section to be taken into account. 1/2
%means the upper alf of the date etc, to remove transient
time_ratio = 1;

%TYPE OF NORMALIZATION -----------------------------

% norma = 1; %If norm 1 then we normalize, if norm 0 then we plot in gauss


%% Choose 1 if wish to plot or zero if not --------------------------------------------
TS = 1; %time series
XC = 0; %cross corr
HG = 0; %histograms
PH = 0; %path
PW = 0; %pwelch
RM = 0; %rms fluctuations
TW = 0; %travelling waves
TM = 0; %thrid moment
SK = 0; %skewness and kurtosis


salva = 0; %Im thinking something to autosave the figures,later
% salvapath='/Users/rubenrojas/Documents/3M/Amplification/3m_smooth/figures_amp/time_series_plots/';
salvapath='/Users/rubenrojas/Documents/3M/Amplification/3m_folder_files/110321';

sz=70; %scatter plots size
fz=20; %font size
lz=20; %legend size
pz=1.5; %plots size


% cd(['/Users/rubenrojas/Documents/3M/Amplification/3m_folder_files/' run(1:6)]);
ampdir = dir([pwd '/' run  'amp.mat']); %seems a directory but its only 1 run
outdir = dir([pwd '/' run  'out.mat']);


load([outdir.folder '/' outdir.name],'folder');

if ~exist('record','var')  %checking if record matrix which is heavy AF is already loaded in workspace
    load([outdir.folder '/' outdir.name],'record');
elseif ~strcmp(record{1,1},run(1:end-2))
    load([outdir.folder '/' outdir.name],'record');
elseif ~strcmp(run,passrun)
    load([outdir.folder '/' outdir.name],'record');
end

load([ampdir.folder '/' ampdir.name],'rm','re','ro1','bext','var');

%  %cleaning spoureous data ro1=-1 or ro1=0
%         a=ro1==0;b=ro1==-1;%c=bext<=1;
%         ro1(a)=NaN; ro1(b)=NaN;%ro1(c)=NaN;
%         bext(a)=NaN; bext(b)=NaN;%bext(c)=NaN;
%         rm(a)=NaN; rm(b)=NaN;%rm(c)=NaN;

if strcmp(path,'full')
    path=1:length(rm);
end



% %prelocting variables
% time=zeros(size(record{k,3}{1,4}(:,32),1),size(path,2));
% brad_t=zeros(size(record{k,3}{1,4}(:,32),1),size(path,2));
% bphi_t=zeros(size(record{k,3}{1,4}(:,32),1),size(path,2));
% coil=zeros(33,size(path,2));

% definning color map____________________________________________________%
cmap=0.8.*[1,0,0.852000000000000;1,0,0;1,0.396000000000000,0;1,0.792000000000000,0;0.812000000000000,1,0;0.416000000000000,1,0;0,1,0.772000000000000;0,0.436000000000000,1;0,0.0400000000000000,1;0.356000000000000,0,1;0.752000000000001,0,1];w=11;
% cmap=[0,0,0;230, 25, 75;245, 130, 48;255, 225, 25;210, 245, 60; 60, 180, 75;70, 240, 240; 0, 130, 200;145, 30, 180;240, 50, 230;250, 190, 212;170, 255, 195;220, 190, 255; 128, 0, 0;170, 110, 40;128, 128, 0; 0, 128, 128;0, 0, 128;128, 128, 128]./255;w=20;
if length(path)>11
cmap=[0,0,0;230, 25, 75;245, 130, 48;255, 225, 25;210, 245, 60; 60, 180, 75;70, 240, 240; 0, 130, 200;145, 30, 180;240, 50, 230;250, 190, 212;170, 255, 195;220, 190, 255; 128, 0, 0;170, 110, 40;128, 128, 0; 0, 128, 128;0, 0, 128;128, 128, 128]./255;w=19;
end
%defining torque funtion for 3m
fun_G_Re_3M = @(l) (l^1.885)*0.003235;

clearvars time brad_t bphi_t coil

%Preeliminares for each point in path------------------------------------------------------------------
for k=path
    %% DEFINNING VARIABLES
    
    if bext(k)<=2   %THIS IS TEMPORALRY WHEN ANALIZYNG B=0 RUNS will change it to normalizing with noise level
        bext(k)=0;
    else
    end
    
    coil=coils_signal_poly(bext(k),var); %this is a 33 vector, being the last ones brad and bphi var is 1 here since it's single coil
    if norma == 1
        if  str2double(folder(5:6)) >= 21
            norm_fac = coil(34);
        else
            norm_fac = coil(32);
        end
    else
        norm_fac = 0.0315;
    end
    %CHANGE time_ratio accordingly to get rid of transient time. Set to 1
    %to take the whole section of data for each parameter of the run.Set to
    %1/2 to choose second half of the parameter ru. etc.
    
    time = record{k,3}{1,2};
    timet = record{k,4}{1,1}; %time for the torque vector, which is not the same as for the Bfield
    time_k = length(record{k,3}{1,2});
    timet_k = length(record{k,4}{1,1});
    time_i = time_k-round(time_k*time_ratio)+1;
    timet_i = timet_k-round(timet_k*time_ratio)+1;
    time = time(time_i:end);
    time = time-min(time)+1; %setting t1=0
    timet =  timet(timet_i:end);
    timet = timet-min(timet)+1; %setting t1=0

    beq1_t = (record{k,3}{1,4}(time_i:end,12)-coil(12))./norm_fac;%radial probes from array
    beq2_t = (record{k,3}{1,4}(time_i:end,13)-coil(13))./norm_fac;%radial probes from array
    beq3_t = (record{k,3}{1,4}(time_i:end,14)-coil(14))./norm_fac;%radial probes from array
    beq4_t = (record{k,3}{1,4}(time_i:end,15)-coil(15))./norm_fac;%radial probes from array
    beq5_t = (record{k,3}{1,4}(time_i:end,16)-coil(16))./norm_fac;%radial probes from array
    beq6_t = (record{k,3}{1,4}(time_i:end,17)-coil(17))./norm_fac;%radial probes from array
    beq7_t = (record{k,3}{1,4}(time_i:end,18)-coil(18))./norm_fac;%radial probes from array
    beq8_t = (record{k,3}{1,4}(time_i:end,19)-coil(19))./norm_fac;%radial probes from array
    beq9_t = (record{k,3}{1,4}(time_i:end,20)-coil(20))./norm_fac;%radial probes from array
    
    bpo1_t = (record{k,3}{1,4}(time_i:end,1)-coil(1))./norm_fac;%radial probes from array
    bpo2_t = (record{k,3}{1,4}(time_i:end,2)-coil(2))./norm_fac;%radial probes from array
    bpo3_t = (record{k,3}{1,4}(time_i:end,3)-coil(3))./norm_fac;%radial probes from array
    bpo4_t = (record{k,3}{1,4}(time_i:end,4)-coil(4))./norm_fac;%radial probes from array

    
    gl1m0_t  =(record{k,3}{1,6}(time_i:end,1));
        gl1mc1_t =(record{k,3}{1,6}(time_i:end,2));
        gl1ms1_t =(record{k,3}{1,6}(time_i:end,3));
    gl1m1_t = (sqrt(record{k,3}{1,6}(time_i:end,2).^2 + record{k,3}{1,6}(time_i:end,3).^2));
    gl2m0_t  =(record{k,3}{1,6}(time_i:end,4));
        gl2mc1_t =(record{k,3}{1,6}(time_i:end,5));
        gl2ms1_t =(record{k,3}{1,6}(time_i:end,6));
    gl2m1_t = (sqrt(record{k,3}{1,6}(time_i:end,5).^2 + record{k,3}{1,6}(time_i:end,6).^2));
        gl2mc2_t =(record{k,3}{1,6}(time_i:end,7));
        gl2ms2_t =(record{k,3}{1,6}(time_i:end,8));
        gl2m2_t = (sqrt(record{k,3}{1,6}(time_i:end,7).^2 + record{k,3}{1,6}(time_i:end,8).^2));
    gl3m0_t  =(record{k,3}{1,6}(time_i:end,9));
        gl3mc1_t =(record{k,3}{1,6}(time_i:end,10));
        gl3ms1_t =(record{k,3}{1,6}(time_i:end,11));
        gl3m1_t = (sqrt(record{k,3}{1,6}(time_i:end,10).^2 + record{k,3}{1,6}(time_i:end,11).^2));
        gl3mc2_t =(record{k,3}{1,6}(time_i:end,12));
        gl3ms2_t =(record{k,3}{1,6}(time_i:end,13));
        gl3m2_t = (sqrt(record{k,3}{1,6}(time_i:end,12).^2 + record{k,3}{1,6}(time_i:end,13).^2));
        gl3mc3_t =(record{k,3}{1,6}(time_i:end,14));
        gl3ms3_t =(record{k,3}{1,6}(time_i:end,15));
        gl3m3_t = (sqrt(record{k,3}{1,6}(time_i:end,14).^2 + record{k,3}{1,6}(time_i:end,15).^2));
    gl4m0_t  =(record{k,3}{1,6}(time_i:end,16));
        gl4mc1_t =(record{k,3}{1,6}(time_i:end,17));
        gl4ms1_t =(record{k,3}{1,6}(time_i:end,18));
        gl4m1_t = (sqrt(record{k,3}{1,6}(time_i:end,17).^2 + record{k,3}{1,6}(time_i:end,18).^2));
        gl4mc2_t =(record{k,3}{1,6}(time_i:end,19));
        gl4ms2_t =(record{k,3}{1,6}(time_i:end,20));
        gl4m2_t = (sqrt(record{k,3}{1,6}(time_i:end,19).^2 + record{k,3}{1,6}(time_i:end,20).^2));
        gl4mc3_t =(record{k,3}{1,6}(time_i:end,21));
        gl4ms3_t =(record{k,3}{1,6}(time_i:end,22));
        gl4m3_t = (sqrt(record{k,3}{1,6}(time_i:end,21).^2 + record{k,3}{1,6}(time_i:end,22).^2));
        gl4mc4_t =(record{k,3}{1,6}(time_i:end,23));
        gl4ms4_t =(record{k,3}{1,6}(time_i:end,24));
        gl4m4_t = (sqrt(record{k,3}{1,6}(time_i:end,23).^2 + record{k,3}{1,6}(time_i:end,24).^2));
    
    %
    %     bg10_t = record{k,3}{1,6}(time_i:end,1).*2./norm_fac;%NOT PROPERLY NORMALIZE Gauss Coefficients here i still need to substract the coils contributions
    %     bg11_t = record{k,3}{1,6}(time_i:end,3).*12./norm_fac;
    %     bg30_t = record{k,3}{1,6}(time_i:end,9).*12./norm_fac;
    %     bg20_t = record{k,3}{1,6}(time_i:end,4).*6./norm_fac;
    %     bg40_t = record{k,3}{1,6}(time_i:end,16).*20./norm_fac;
    
    if norma == 1
            
            
            if bext(k)<1  %this is in case of 0 external applied field
                coil=coils_signal_poly(20,var);
                fprintf('Element path = %1.0f normalized with respect to 20 amps in the coils.\n',i)

            else
                coil=coils_signal_poly(bext(i),var); %this is a 33 vector, being the last ones brad and bphi
            end
                coil_coeff=gcoeff3m(coil(1:31),probepos);%this is a 24 vector being the 24 gauss coeff l=4
%                 coil_coeff = zeros(1,24); %doing this if gainplot2 was used
                
                norm_fac = coil(32);
          
            
        elseif norma == 0
%           fprintf('Bfield in Gauss.\n')
            norm_fac = 0.0315;
            if bext(k)<1
            coil_coeff = zeros(1,24); %this was first attemp of plotting without substracting coils
            coil = zeros(1,35);
            else
        
            coil=coils_signal_poly(bext(k),var); %this is a 33 vector, being the last ones brad and bphi
%              coil_coeff = zeros(1,24); %doing this if gainplot2 was used
            coil_coeff=gcoeff3m(coil(1:31),probepos);%this is a 24 vector being the 24 gauss coeff l=4

            end
            
        
    end
        
        if str2double(folder(5:6)) >= 21
            fun_G_Re_3M = @(l) (l^1.845416072636967)*0.041601147102921+1.313735481971338e+11;
        else
            fun_G_Re_3M = @(l) (l^1.885)*0.003235;
        end
      
        
      
        %compensating for choosing different probe for runs between nov2021
        %and feb 2022____________________________________________________%
        if str2double(folder(5:6)) == 21 || str2double(folder(5:6)) == 22
            if str2double(folder(1:2)) == 11 || str2double(folder(1:2)) == 12 || str2double(folder(1:2)) == 01|| str2double(folder(1:2)) == 02
              if norma == 0
                    norm_fac_phi = 0.005; %sensitivity of the phi(-) probe
              elseif norma == 1
                    norm_fac_phi = norm_fac/6.3; %ratio of sensitivity of the different probes
              end
            else
               norm_fac_phi = norm_fac;
            end
        else
          norm_fac_phi = norm_fac;
        end
    
    brad_t = (record{k,3}{1,4}(time_i:end,32)-coil(32))./norm_fac;
    bphi_t = (record{k,3}{1,4}(time_i:end,33)-coil(33))./norm_fac_phi;
    if str2double(folder(5:6)) >= 21
    brad2_t = (record{k,3}{1,4}(time_i:end,34)-coil(34))./norm_fac;
    bphi2_t = (record{k,3}{1,4}(time_i:end,35)-coil(35))./norm_fac;
    end
    
    bl1m0_t  =(2*0.95^3*(gl1m0_t-coil_coeff(1))/norm_fac);
        bl1mc1_t =2*0.95^3*(gl1mc1_t-coil_coeff(2))/norm_fac;
        bl1ms1_t =2*0.95^3*(gl1ms1_t-coil_coeff(3))/norm_fac;
    bl1m1_t= 2*0.95^3*(gl1m1_t-sqrt(coil_coeff(2)^2+coil_coeff(3)^2))/norm_fac;
    bl2m0_t  =(6*0.95^4*(gl2m0_t-coil_coeff(4))/norm_fac);
        bl2mc1_t =6*0.95^4*(gl2mc1_t-coil_coeff(5))/norm_fac;
        bl2ms1_t =6*0.95^4*(gl2ms1_t-coil_coeff(6))/norm_fac;
    bl2m1_t=6*0.95^4*(gl2m1_t-sqrt(coil_coeff(5)^2+coil_coeff(6)^2))/norm_fac;
        bl2mc2_t =6*0.95^4*(gl2mc2_t-coil_coeff(7))/norm_fac;
        bl2ms2_t =6*0.95^4*(gl2ms2_t-coil_coeff(8))/norm_fac;
        bl2m2_t=6*0.95^4*(gl2m2_t-sqrt(coil_coeff(7)^2+coil_coeff(8)^2))/norm_fac;
    bl3m0_t  =(12*0.95^5*(gl3m0_t-coil_coeff(9))/norm_fac);
        bl3mc1_t =12*0.95^5*(gl3mc1_t-coil_coeff(10))/norm_fac;
        bl3ms1_t =12*0.95^5*(gl3ms1_t-coil_coeff(11))/norm_fac;
        bl3m1_t=12*0.95^5*(gl3m1_t-sqrt(coil_coeff(10)^2+coil_coeff(11)^2))/norm_fac;
        bl3mc2_t =12*0.95^5*(gl3mc2_t-coil_coeff(12))/norm_fac;
        bl3ms2_t =12*0.95^5*(gl3ms2_t-coil_coeff(13))/norm_fac;
        bl3m2_t=12*0.95^5*(gl3m2_t-sqrt(coil_coeff(13)^2+coil_coeff(12)^2))/norm_fac;
        bl3mc3_t =12*0.95^5*(gl3mc3_t-coil_coeff(14))/norm_fac;
        bl3ms3_t =12*0.95^5*(gl3ms3_t-coil_coeff(15))/norm_fac;
        bl3m3_t=12*0.95^5*(gl3m3_t-sqrt(coil_coeff(14)^2+coil_coeff(15)^2))/norm_fac;
    bl4m0_t  =(20*0.95^6*(gl4m0_t-coil_coeff(16))/norm_fac);
        bl4mc1_t =20*0.95^6*(gl4mc1_t-coil_coeff(17))/norm_fac;
        bl4ms1_t =20*0.95^6*(gl4ms1_t-coil_coeff(18))/norm_fac;
        bl4m1_t = 20*0.95^6*(gl4m1_t-sqrt(coil_coeff(17)^2+coil_coeff(18)^2))/norm_fac;
        bl4mc2_t =20*0.95^6*(gl4mc2_t-coil_coeff(19))/norm_fac;
        bl4ms2_t =20*0.95^6*(gl4ms2_t-coil_coeff(20))/norm_fac;
        bl4m2_t=20*0.95^6*(gl4m2_t-sqrt(coil_coeff(19)^2+coil_coeff(20)^2))/norm_fac;
        bl4mc3_t =20*0.95^6*(gl4mc3_t-coil_coeff(21))/norm_fac;
        bl4ms3_t =20*0.95^6*(gl4ms3_t-coil_coeff(22))/norm_fac;
        bl4m3_t=20*0.95^6*(gl4m3_t-sqrt(coil_coeff(21)^2+coil_coeff(22)^2))/norm_fac;
        bl4mc4_t =20*0.95^6*(gl4mc4_t-coil_coeff(23))/norm_fac;
        bl4ms4_t =20*0.95^6*(gl4ms4_t-coil_coeff(24))/norm_fac;
        bl4m4_t=20*0.95^6*(gl4m4_t-sqrt(coil_coeff(23)^2+coil_coeff(24)^2))/norm_fac;
    %
    p1_t = record{k,3}{1,3}(time_i:end,1)./4;
    p2_t = record{k,3}{1,3}(time_i:end,2)./4;
    p3_t = record{k,3}{1,3}(time_i:end,3)./4;
    
    %Torque is not at 256 so it cannot be plotted side to side with mag
    %data, it has its own plotter
    torinn_t=record{k,4}{1,2}(:,1);
    g_t=abs(torinn_t./(((0.7e-6)^2)*(0.45)*927).*(1130/1801990));%constant found in Doug codes used to normalized torque data. whERE THE heck comes from, IDK
    gginf_t=g_t./fun_G_Re_3M(re(k));
    gginf_t =  gginf_t(time_i:end);
    g_t = g_t(time_i:end);
    time_torque=record{k,4}{1,1}(:,1);
    time_torque= time_torque(time_i:end);
    
    %Normalizing Time from sec to outer_rot_rate. NOTE:Since we are setting
    %time t1=0 the PHASE IS UNKNOWN
    
    out_freq=abs(record{k,2}(1,2)); %outer frequency in Hz for normalice
    if out_freq <= 0.001
        out_freq = 1;
    fprintf('Outer sphere freq set to unity\n')

    end
    inn_freq=abs(record{k,2}(1,3)); %inner frequency for normalice
    
%     
    time = time*out_freq;
    timet = timet*out_freq;
    
    %% TIME SERIES
    if TS ==1
        
        max_lag=11; %in rev of the outer
        
        figure(45);
        if k==path(1);clf;end
        
        subplot(2,1,1)
        plot(timet,movmean(gginf_t,30*1),'Color',cmap(mod(k+3,10)+1,:),'LineWidth',pz);hold on;
%         plot(timet,gginf_t,'Color',cmap(mod(k+3,10)+1,:),'LineWidth',pz);hold on;
            
            xlabel('$\Omega_{o}t/2 \pi$','FontName','Times','Interpreter','latex');
            ylabel('$G/G_{\infty} $','FontName','Times','Interpreter','latex');            
            if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            set(gca,'linewidth',2,'FontName','Times New Roman','FontSize',fz);
            hold off
% xlim([1000 1500])
title([run(1:6) ' | Path=' num2str(path)],'interpreter','latex');

%         figure(46);
        subplot(2,1,2)
%         plot(time,bphi_t,'Color',cmap(mod(k+3,10)+1,:),'LineWidth',pz);hold on;
        plot(time,(movmean(bl1m0_t,256*1)),'color',cmap(mod(k+3,10)+1,:),'LineWidth',pz);hold on;
%         xlim([900 1600])
        %         plot(time,brad2_t,'color',cmap(mod(k+3,10)+1,:),'LineWidth',1);
%         plot(time,bphi2_t,'color',cmap(mod(k+7,10)+1,:),'LineWidth',1);
        
        %                         plot(time,beq1_t,'Color',cmap(mod(k,w)+1,:),'LineWidth',1);
        
        if k==path(end)
           
            %             xlim([(max(time)-min(time))/2 min(time)+(max(time)-min(time))/2+5])
            %             xlim([10 max_lag]) %setting up to
                        xlabel('$\Omega_{o}t/2 \pi$','FontName','Times','Interpreter','latex');
%             xlabel('$t (\mathrm{ssm}) $','FontName','Times','Interpreter','latex');
                        ylabel('$B_{\phi}$','FontName','Times','Interpreter','latex');
%             ylabel('$B_{\phi}[G]$','FontName','Times','Interpreter','latex');
            
            if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            set(gca,'linewidth',2,'FontName','Times New Roman','FontSize',fz);
            hold off
%             xlim([1000 1500])
            
            if salva==1
                saveas(gcf,[salvapath '/time_series_brad.png'],'png');
            else
            end
            
        end
        %
        %         figure(46); hold on;
        %         if k==path(1);clf;end
        %         %         plot(time,brad_t,'Color',cmap(mod(k,w)+1,:),'LineWidth',1);
        %         plot(time,bphi_t,'color',cmap(mod(k,w)+1,:),'LineWidth',1);
        %         if k==path(end)
        % %             xlim([(max(time)-min(time))/2 min(time)+(max(time)-min(time))/2+5])
        %             xlim([1 max_lag]) %setting up to 3 rev of the outer
        %             xlabel('$\Omega_{o}t/2 \pi$','FontName','Times','Interpreter','latex');
        %             ylabel('$B_{\varphi}/B_0$','FontName','Times','Interpreter','latex');
        %             if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
        %             if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
        %             if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
        %             set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
        %             hold off
        %
        %
        %         if salva==1
        %             saveas(gcf,[salvapath '/time_series_bphi.png'],'png');
        %         else
        %         end
        %
        %         end
        %
        %         figure(47); hold on;
        %         if k==path(1);clf;end
        %         plot(time,beq1_t,'Color',cmap(mod(k,w)+1,:),'LineWidth',1);
        %         if k==path(end)
        % %             xlim([(max(time)-min(time))/2 min(time)+(max(time)-min(time))/2+5])
        %             xlim([1 max_lag]) %setting up to 3 rev of the outer
        %             xlabel('$\Omega_{o}t/2 \pi$','FontName','Times','Interpreter','latex');
        %             ylabel('$B_{eq1}/B_0$','FontName','Times','Interpreter','latex');
        %             if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
        %             if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
        %             if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
        %             set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
        %             hold off
        %
        %
        %         if salva==1
        %             saveas(gcf,[salvapath '/time_series_beq1.png'],'png');
        %         else
        %         end
        %
        %         end
        
        
        
        
    else
    end
    
    %% XCORR
    if XC ==1
        figure(75);hold on
        if k==path(1);clf;end
        
        %select maximum lag in sec
        %         max_lag=10; %in secs
        max_lag = 10000; %in rev of outer sphere
        %calculating xcorr for finger: HELP xcorr for more options
        %'normalized' or 'coeff' ? Normalizes the sequence so that the
        %autocorrelations at zero lag equal 1:
%         [r,lags]=xcorr(brad_t,bphi_t,round(max_lag*256./out_freq),'normalized');
%         [r,lags]=xcorr(brad_t-mean(brad_t),bphi_t-mean(bphi_t),round(max_lag*256./inn_freq),'normalized');

        [r,lags]=xcorr(movmean(bphi_t-mean(bphi_t), 256*10),interp1(timet,movmean(g_t-mean(g_t),30*10),time),round(max_lag*30./out_freq),'normalized');
        
        
        %normalicing techniche
        r_nor=1;
%                 r_nor=xcorr(brad_t,brad_t,round(max_lag*256./out_freq),'normalized');
%                 r_nor = r(lags==0);
        
        %plotting
        plot((lags./256).*inn_freq,r./r_nor,'Color',cmap(mod(k,w)+1,:),'LineWidth',2);
        
        if k==path(end)
            %             xlim([-10 10]); %zooming region of the xcorr
            xlabel('$\Omega_{o} \Delta t/2 \pi$','FontName','Times','Interpreter','latex');
            set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
            if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            %             line([0 0], [0.4 1],'Color','k','LineStyle','--','LineWidth',1,'handlevisibility','off');
            %                         ylim([0.4 1]);
            title('$Xcorr(B_{rad},B_{phi})$','interpreter','latex');
            box on
            hold off
            
            
            if salva==1
                saveas(gcf,[salvapath '/xcorre.png'],'png');
            else
            end
            
        end
        
    else
    end
    
    %% HISTOGRAM
    if HG==1
        
         clearvars hall_t
        hall_t = [brad2_t'; bphi2_t'; brad_t';bphi_t']';
%         hall_t = [bl1m0_t';bl3m0_t';bl2m0_t';bl1m1_t']';

        
        %__________________________________________________________________%
        %loop for plotting
        
        
        for i = 1:size(hall_t,2)
            figure(79+i);hold on
            if k==path(1);clf;end
     
        %check binmetohd and normalization in function help:Scott?s rule is
        %optimal if the data is close to being normally distributed. This
        %rule is appropriate for most other distributions, as well. It uses
        %a bin width of 3.5*std(X(:))*numel(X)^(-1/3).
        
        %This line is and example of for centered histograms
                histogram((hall_t(:,i)-mean(hall_t(:,i)))./rms(hall_t(:,i)),'binmethod','scott','DisplayStyle','stairs','normalization','pdf','EdgeColor',cmap(mod(k,w)+1,:),'LineWidth',pz);
%         
%         histogram(hall_t(:,i),'binmethod','scott','DisplayStyle','stairs','normalization','pdf','EdgeColor',cmap(mod(k,w)+1,:),'LineWidth',2);
        
        
        if k==path(end)
            ylabel('${\rm PDF}$','FontName','Times','Interpreter','latex');
            %             xlabel('$B_{r}/B_0$','FontName','Times','Interpreter','latex');
                        xlabel('$(B-\left < B_r \right >) /B_{\rm rms}$','FontName','Times','Interpreter','latex');
%             xlabel(['$B_{' num2str(i) '}(\mathrm{G})$'],'FontName','Times','Interpreter','latex');
           
            if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            set(gca,'linewidth',2,'FontName','Times New Roman','FontSize',fz,'yscale','log');
          title(['$B_{' num2str(i) '}$'],'interpreter','latex');
            ylim([0.001 10]);

            box on
            hold off
            
            if salva==1
                saveas(gcf,[salvapath '/histo_brad.png'],'png');
            else
            end
            
        end
%         
%         figure(77);hold on
%         if k==path(1);clf;end
%         %This line is and example of for centered histograms
%         %         histogram((bphi_t-mean(bphi_t))./rms(bphi_t),'binmethod','scott','DisplayStyle','stairs','normalization','pdf','EdgeColor',cmap(mod(k,w)+1,:),'LineWidth',pz);
%         
%         histogram(bphi2_t,'binmethod','scott','DisplayStyle','stairs','normalization','pdf','EdgeColor',cmap(mod(k,w)+1,:),'LineWidth',2);
%         
%         if k==path(end)
%             ylabel('${\rm PDF}$','FontName','Times','Interpreter','latex');
%             %             xlabel('$(B_{\varphi}-\left < B_{\varphi} \right >) /B_{\rm rms}$','FontName','Times','Interpreter','latex');
%             %             xlabel('$B_{\varphi}/B_0$','FontName','Times','Interpreter','latex');
%             xlabel('$B_{\varphi}(\mathrm{G})$','FontName','Times','Interpreter','latex');
%             
%             set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
%             if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
%             if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
%             if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
%             set(gca,'linewidth',2,'FontName','Times New Roman','FontSize',fz,'yscale','log');
%             box on
%             hold off
%             
%             if salva==1
%                 saveas(gcf,[salvapath '/histo_bphi.png'],'png');
%             else
%             end
%             
%         end
%         
%         figure(78);hold on
%         if k==path(1);clf;end
%         
%         %This line is and example of for centered histograms
%         %         histogram((beq1_t-mean(beq1_t))./rms(beq1_t),'binmethod','scott','DisplayStyle','stairs','normalization','pdf','EdgeColor',cmap(mod(k,w)+1,:),'LineWidth',pz);
%         
%         histogram(beq1_t,'binmethod','scott','DisplayStyle','stairs','normalization','pdf','EdgeColor',cmap(mod(k,w)+1,:),'LineWidth',2);
%         
%         if k==path(end)
%             ylabel('${\rm PDF}$','FontName','Times','Interpreter','latex');
%             %             xlabel('$(B_{eq}-\left < B_{eq} \right >) /B_{\rm rms}$','FontName','Times','Interpreter','latex');
%             %             xlabel('$B_{eq}/B_0$','FontName','Times','Interpreter','latex');
%             xlabel('$B_{eq}(\mathrm{G})$','FontName','Times','Interpreter','latex');
%             
%             set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
%             if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
%             if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
%             if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
%             set(gca,'linewidth',2,'FontName','Times New Roman','FontSize',fz,'yscale','log');
%             box on
%             hold off
%             
%             if salva==1
%                 saveas(gcf,[salvapath '/histo_beq1.png'],'png');
%             else
%             end
        
        
    
        end
    
    end
    
    %% PATH
    if PH==1
        t11=50000;
        t22=60000;
        td=100;
        figure(18);
        plot(brad_t(t11:t22),bphi_t(t11+td:t22+td),'linewidth',1);
        %         scatter(brad_t(t11:t22),bphi_t(t11+td:t22+td),2);
        xlabel('$B_r(t)$','FontName','Times','Interpreter','latex');
        ylabel('$B_{r}(t+1)$','FontName','Times','Interpreter','latex');
        box on
        set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
        hold off
        if salva==1
            saveas(gcf,['/Users/rubenrojas/Documents/3M/Amplification/3m_smooth/figures_amp/time_series_plots/pathrad'],'png');
        else
        end
        %
        % figure(19);
        % % plot(bphi_t(t11:t22),bphi_t(t11+td:t22+td),'linewidth',1);
        % scatter(bphi_t(t11:t22),bphi_t(t11+td:t22+td),1);
        % xlabel('$B_{\varphi}(t)$','FontName','Times','Interpreter','latex');
        % ylabel('$B_{\varphi}(t+1)$','FontName','Times','Interpreter','latex');
        % box on
        % set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
        % hold off
        % if salva==1
        %     saveas(gcf,['/Users/rubenrojas/Documents/3M/Amplification/3m_smooth/figures_amp/time_series_plots/pathphi'],'png');
        % else
        % end
        
    else
    end
  %% PWELCH
    if PW==1
        
        if k==path(1);p_rad=zeros(size(path,2),1);end
        out_freq=abs(record{k,2}(1,2)); %outer frequency for normalice
        if out_freq <= 0.001
            out_freq = 1;
        end
        inn_freq=abs(record{k,2}(1,3)); %inner frequency for normalice
        
        clearvars ball_t
%         ball_t = [brad_t'; bphi_t' ;beq1_t';bl4m3_t';bl1m1_t']';
%          ball_t = detrend(lowpass(brad2_t,60,256));
           ball_t = (brad_t);
        %__________________________________________________________________%
        %loop for plotting
        

        for i = 1:size(ball_t,2)
            figure(101+i);
            if k==path(1);cla;a=0;end;hold on
            
            set(gcf,'color','w','position',[136    46   438   651]);

           
            %__________________________________________________________________%
            %alternatively with default values
%             [pxx,f] = pwelch(ball_t(:,i));
            
            
            
            %%by hand pwelch parameters
%             [pxx,f] = pwelch(ball_t(:,i),16385,16384/2,16385,256);
              [pxx,f] = pwelch(ball_t(:,i),[],[],[],256);
           
            
            %__________________________________________________________________%
            %plotting
            
%              plot(f./out_freq,10*log10(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);%indb
%              plot((f),(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);hold on
%              plot((f./inn_freq),(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',pz);hold on
              
               a=1+a;
%                plot((f./out_freq),(pxx).*10^(k-1),'Color',cmap(a+1,:),'LineWidth',pz);hold on
               plot((f./out_freq),(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',pz);hold on


            
            % Divide by this for compensated plots ./(f./out_freq).^(-11/3)
            
            
            
            %         xlim([0 50]);
            %         axis([10^-2 10^2 10^-5 10^2]);
            %         ylim([-40 20]);
            %         line([out_freq out_freq], [-100 100],'Color','k','LineStyle','--','LineWidth',1,'handlevisibility','off');
            %         line([inn_freq inn_freq], [-100 100],'Color','k','LineStyle','--','LineWidth',1,'handlevisibility','off');
            
            %__________________________________________________________________%
            %calculate slope of spectrum from within given range
            im = 3; % in Out Freq units (see plot for determine)
            fm = 40;
            j=f>im*out_freq & f<fm*out_freq; %tranform from out rot rate to hertz cuz the slope is calculate in hz
            
            [pk,S] = polyfit(log10(f(j)), log10(pxx(j)),1);
            p_rad(k,1)=pk(1);
            p_rad(k,2)=1 - (S.normr/norm(log10(pxx(j)) - mean(log10(pxx(j)))))^2;
            
            %                 p = x.^(-11/3);
            %                 p = polyfit(f(j), pxx(j),1);
            %                 p = polyfit(f(j), pxx(j),1);
            %                 line([(im) (fm)], [(p(1)*log10(im)+p(2)) (p(1)*log10(fm)+p(2))],'Color','k','LineStyle','--','LineWidth',2) % for decibelts
            
            
            %                 p = polyfit(x,y,1)
            %                 yfit = polyval(pk,log10(f(j)));
            %                 yfit =  pk(1) * log10(f(j)) + pk(2);
            %                 yresid = log10(pxx(j)) - yfit;
            %                 SSresid = sum(yresid.^2);
            %                 SStotal = (length(log10(pxx(j)))-1) * var(log10(pxx(j)));
            %                 rsq = 1 - SSresid/SStotal
            
            %__________________________________________________________________%
            %legend and make up for the graph
            if k==path(end)
                %             m=plot((im:fm),(im:fm).^(-11/3)*10^(pk(2)),'Color','k','LineStyle','--','LineWidth',2,'handlevisibility','off');
                xlabel('$\omega/\Omega_o$','FontName','Times','Interpreter','latex');
                %             ylabel('dB','FontName','Times','Interpreter','latex');
                if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','southwest','Interpreter','latex');else;end
                if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','southwest','Interpreter','latex');else;end
                if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','southwest','Interpreter','latex');else;end
                set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz,'xscale','log','yscale','log');
                box on
                
%                 title(['$B_{' num2str(i) '}$'],'interpreter','latex');
                
                hold off
                
                
                if salva==1
                    saveas(gcf,[salvapath '/spec_brad.png'],'png');
                else
                end
                
            end %end on legend labels etc
            
        end %end of plotting
        %         %         % SECOND FIGURES
        %
        %         figure(94);hold on
        %         if k==path(1);clf;p_phi=zeros(size(path,2),1);end
        %         out_freq=abs(record{k,2}(1,2)); %outer frequency for normalice
        %         inn_freq=abs(record{k,2}(1,3)); %inner frequency for normalice
        %
        %         %alternatively with default values
        %         %         [pxx,f] = pwelch(brad_t);
        %
        %
        %         %%by hand pwelch parameters
        %         %         [pxx,f] = pwelch(bphi2_t,16385,16384/2,16385,256);
        %         [pxx,f] = pwelch(bphi2_t);
        %         %                 plot(f./out_freq,10*log10(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);%in db
        %         %         plot(f,(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);%in db
        %         plot((f),(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);hold on
        %
        %         %         plot(f./out_freq,(pxx)./(f./out_freq).^(-11/3),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);
        %
        %         %         xlim([0 50]);
        %         %         axis([10^-2 10^2 10^-5 10^2]);
        %         %         ylim([-40 20]);
        %         %         line([out_freq out_freq], [-100 100],'Color','k','LineStyle','--','LineWidth',1,'handlevisibility','off');
        %         %         line([inn_freq inn_freq], [-100 100],'Color','k','LineStyle','--','LineWidth',1,'handlevisibility','off');
        %
        %         %        %calculate slope of spectrum from within given range
        %         im = 4; % in Out Freq units (see plot for determine)
        %         fm = 20;
        %         j=f>im*out_freq & f<fm*out_freq; %tranform from out rot rate to hertz cuz the slope is calculate in hz
        %
        %
        %         [pk,S] = polyfit(log10(f(j)), log10(pxx(j)),1);
        %         p_phi(k,1)=pk(1);
        %         p_phi(k,2)=1 - (S.normr/norm(log10(pxx(j)) - mean(log10(pxx(j)))))^2;
        %
        %         %                 p = polyfit(f(j), pxx(j),1);
        %         %                 line([(2) (20)], [(p(1)*log10(2)+p(2)) (p(1)*log10(20)+p(2))],'Color','k','LineStyle','--','LineWidth',2) % for decibelts
        %
        %
        %         if k==path(end)
        %             xlabel('$\omega/\Omega{o}$','FontName','Times','Interpreter','latex');
        %             %             ylabel('dB','FontName','Times','Interpreter','latex');
        %             if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','southwest','Interpreter','latex');else;end
        %             if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','southwest','Interpreter','latex');else;end
        %             if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','southwest','Interpreter','latex');else;end
        %             set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz,'xscale','log','yscale','log');
        %             box on
        %             title('$B_{\varphi}$','interpreter','latex');
        %
        %             hold off
        %
        %
        %             if salva==1
        %                 saveas(gcf,[salvapath '/spec_bphi.png'],'png');
        %             else
        %             end
        %         end
        %
        %
        %         figure(95);hold on
        %         if k==path(1);clf;p_eq=zeros(size(path,2),1);end
        %         out_freq=abs(record{k,2}(1,2)); %outer frequency for normalice
        %         inn_freq=abs(record{k,2}(1,3)); %inner frequency for normalice
        %
        %         %alternatively with default values
        % %         [pxx,f] = pwelch(beq1_t);
        %
        %
        %         %         %%by hand pwelch parameters
        %                 [pxx,f] = pwelch(beq1_t,16385,16384/2,16385,256);
        %
        %         %                 plot(f./out_freq,10*log10(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);%indb
        % %                 plot(f./out_freq,(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);%indb
        %         plot((f./inn_freq),(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);hold on
        % %                         plot(f,(pxx),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);%indb
        %
        %         %         plot(f./out_freq,(pxx)./(f./out_freq).^(-11/3),'Color',cmap(mod(k,w)+1,:),'LineWidth',1);
        %
        %         %         xlim([0 50]);
        %         %         ylim([-40 20]);
        %         %         line([out_freq out_freq], [-100 100],'Color','k','LineStyle','--','LineWidth',1,'handlevisibility','off');
        %         %         line([inn_freq inn_freq], [-100 100],'Color','k','LineStyle','--','LineWidth',1,'handlevisibility','off');
        %
        %         %         %calculate slope of spectrum from within given range
        %         im = 1; % in Out Freq units (see plot for determine)
        %         fm = 2;
        %         j=f>im*out_freq & f<fm*out_freq; %tranform from out rot rate to hertz cuz the slope is calculate in hz
        %
        %         [pk,S] = polyfit(log10(f(j)), log10(pxx(j)),1);
        %         p_eq(k,1)=pk(1);
        %         p_eq(k,2)=1 - (S.normr/norm(log10(pxx(j)) - mean(log10(pxx(j)))))^2;
        %
        %         %                 p = polyfit(f(j), pxx(j),1);
        %         %                 line([(2) (20)], [(p(1)*log10(2)+p(2)) (p(1)*log10(20)+p(2))],'Color','k','LineStyle','--','LineWidth',2) % for decibelts
        %
        %
        %         if k==path(end)
        %             xlabel('$\omega/\Omega{o}$','FontName','Times','Interpreter','latex');
        %             %             ylabel('dB','FontName','Times','Interpreter','latex');
        %             if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','southwest','Interpreter','latex');else;end
        %             if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','southwest','Interpreter','latex');else;end
        %             if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','southwest','Interpreter','latex');else;end
        %             set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz,'xscale','log','yscale','log');
        %             box on
        %             title('$B_{eq}$','interpreter','latex');
        %
        %             hold off
        %
        %
        %             if salva==1
        %                 saveas(gcf,[salvapath '/spec_beq.png'],'png');
        %             else
        %             end
        %
        %         end
        
        
        
        
    end %end of section
    %% RMS fluctuaciones
    if RM==1
        figure(99);hold on
        if k==path(1);clf;end
        
        %Here you choose what to plot
        
%        ball_t = [beq1_t'; beq2_t' ;beq3_t';beq4_t';beq5_t';beq6_t'; beq7_t' ;beq8_t']';
       ball_t = [brad_t';brad2_t';bphi2_t']';
%        ball_t = [bpo1_t'; bpo2_t' ;bpo3_t';bpo4_t']';

        for i = 1:size(ball_t,2)
           b=ball_t(:,i);
        
        if strcmp(ramp,'ros');scatter(ro1(k),std(b),sz,'filled','markerfacecolor',cmap(i,:));xlabel('$Ro^{-1}$','FontName','Times','Interpreter','latex');hold on;end
        if strcmp(ramp,'rey');scatter(rm(k),std(b),sz,'filled','markerfacecolor',cmap(i,:));xlabel('$Rm$','FontName','Times','Interpreter','latex');hold on;end
        if strcmp(ramp,'mag');scatter(bext(k),std(b),sz,'k','filled');xlabel('$B_{ext}$','FontName','Times','Interpreter','latex');hold on;end
        
         
        
        
        if k==path(end)
            ylabel('$\sigma_{B_{eq}}/<|B_{eq}|>$','FontName','Times','Interpreter','latex');
            set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
%             title('$B_{eq}$', 'interpreter','latex');
            %          ylim([0.2 1.8]);
            box on
%             hold off
            
            
            if salva==1
                saveas(gcf,[salvapath '/rms_eq.png'],'png');
            else
            end
            
        end
        end
    end
    
    
    
    %% TRAVELING WAVES
    if TW ==1
        inni=100; %initial point
        max_lag=1; %in rev of the outer
        
        figure(44);
        if k==path(1);clf;end
        
        
        tw1=plot(time,beq1_t,'Color',cmap(mod(k+1,w),:),'LineWidth',2);hold on;
        tw2=plot(time,beq2_t,'color',cmap(mod(k+1,w),:),'LineWidth',2);
        tw3=plot(time,beq3_t,'color',cmap(mod(k+3,w),:),'LineWidth',2);
        
        %                         plot(time,beq1_t,'Color',cmap(mod(k,w)+1,:),'LineWidth',1);
        
        if k==path(end)
            %             xlim([(max(time)-min(time))/2 min(time)+(max(time)-min(time))/2+5])
            xlim([inni inni+max_lag]) %setting up to
            xlabel('$\Omega_{o}t/2 \pi$','FontName','Times','Interpreter','latex');
            ylabel('$B_{r}/B_0$','FontName','Times','Interpreter','latex');
            if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
            hold off
            
            legend([tw1,tw2,tw3],{'$1$','$2$','$3$'},'location','northeast','Interpreter','latex');
            
            
            
            if salva==1
                saveas(gcf,[salvapath '/time_series_brad.png'],'png');
            else
            end
        end
    end
    
    
    %% Thrid MOMMENT
    if TM == 1
        figure(24);hold on
        if k==path(1);clf;m_third=zeros(size(path,2),1);end
        
        out_freq=abs(record{k,2}(1,2)); %outer frequency for normalice
        %         inn_freq=abs(record{k,2}(1,3)); %inner frequency for normalice
        
        lag = 0.3; %in units of outer rot or inner CHECK Below
        lags = round(lag*256./out_freq);
        %         lags = round(lag*256./inn_freq);
        
        vec = brad_t(1:end); %choosing vector to do the thing to. also choosign also a pice of 10 times the size of the lag to multiple
        lag_t = (1:lags);
        r = zeros(1,length(lag_t));
        i=0;
        for t = 1:lags
            i = 1+i;
            r(i) = mean((vec(1:end-t) - vec(t+1:end)).^3);
        end
        %
        %         r = [1:10000];
        %         r2 = r.^2;
        
        %plotting
        plot((lag_t./256).*out_freq,r,'Color',cmap(mod(k,w)+1,:),'LineWidth',2);
        %         plot(r,r2);
        
        %calculate slope of spectrum from within given range
        if ismember(k,16:24)
            im = 0.01; % in Out Freq units (see plot for determine)
            fm = 0.1;
        else
            im = 0.01; % in Out Freq units (see plot for determine)
            fm = 0.04;
        end
        j=(lag_t./256).*out_freq>im*out_freq & (lag_t./256).*out_freq<fm*out_freq; %tranform from out rot rate to hertz cuz the slope is calculate in hz
        %                 j = r>im & r<fm;
        %                 pk = polyfit(log10(r), log10(r2),1);
        
        pk = polyfit(log10((lag_t(j)./256).*out_freq), log10(r(j)),1);
        m_third(k)=pk(1);
        
        %                 line([(im) (fm)], [(p(1)*(im)+p(2)) (p(1)*(fm)+p(2))],'Color','k','LineStyle','--','LineWidth',2,'handlevisibility','off') % for decibelts
        if k == path(end)
            %                     line([im fm], [0.3 (10^3)*(fm-im)+0.3],'Color','k','LineStyle','--','LineWidth',2,'handlevisibility','off')
            %                     line([im fm], [0.7 (10^3)*(fm-im)+0.7],'Color','k','LineStyle','--','LineWidth',2,'handlevisibility','off')
            %                 p = polyfit((lag_t(j)./256).*out_freq, r(j),1);
            %                 line([(im) (fm)], [(p(1)*(im)+p(2)) (p(1)*(fm)+p(2))],'Color','k','LineStyle','--','LineWidth',2,'handlevisibility','off') % for decibelts
        end
        
        
        
        if k==path(end)
            %             xlim([-10 10]); %zooming region of the xcorr
            xlabel('$\Omega_{o} s/2 \pi$','FontName','Times','Interpreter','latex');
            set(gca,'linewidth',2,'FontName','Times New Roman','FontSize',fz);
            set(gca,'xscale', 'log', 'yscale','log');
            if strcmp(ramp,'ros');leg=cell(1,length(path));leg(:)={'$Ro^{-1}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(ro1(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'rey');leg=cell(1,length(path));leg(:)={'$Rm='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(rm(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            if strcmp(ramp,'mag');leg=cell(1,length(path));leg(:)={'$B_{ext}='};leg2=cell(1,length(path));leg2(:)={'$'};lgd=legend(strcat(leg',num2str(bext(path)',3),leg2'),'fontsize',lz,'location','northwest','Interpreter','latex');else;end
            %             line([0.002 0.03], [(-2*log10(0.002)) (-2*log10(0.03))],'Color','k','LineStyle','-','LineWidth',2,'handlevisibility','off');
            %             line([0.002 0.03], [2*0.002+0.7 2*0.03+0.7],'Color','k','LineStyle','-','LineWidth',2,'handlevisibility','off');
            %                         ylim([0.4 1]);
            title('$\left < \left(B(t+s) - B(t) \right) ^3 \right >$','interpreter','latex');
            box on
            hold off
            
            
            if salva==1
                saveas(gcf,[salvapath '/xcorre.png'],'png');
            else
            end
            
        end
        
    else
    end
    
    
    %% Skewness and Kurtosis
    if SK==1
        figure(46);hold on
        if k==path(1);clf;end
        
        %Here you choose what to plot
        b=bphi_t;
        
        if strcmp(ramp,'ros');scatter(ro1(k),skewness(b),sz,'k','filled');xlabel('$Ro^{-1}$','FontName','Times','Interpreter','latex');hold on;end
        if strcmp(ramp,'rey');scatter(rm(k),skewness(b),sz,'k','filled');xlabel('$Rm$','FontName','Times','Interpreter','latex');hold on;end
        if strcmp(ramp,'mag');scatter(bext(k),skewness(b),sz,'k','filled');xlabel('$B_{ext}$','FontName','Times','Interpreter','latex');hold on;end
        
        
        
        if k==path(end)
            ylabel('$Sk$','FontName','Times','Interpreter','latex');
            set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
            title('$B_{}$', 'interpreter','latex');
            %          ylim([0.2 1.8]);
            box on
            hold off
        end
    end
    
    
end %the whole thing

