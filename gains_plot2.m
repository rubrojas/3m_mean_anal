%gains_plot2 Ruben Feb 2022
%uses record matrix from Artur
%This is to work with a single record matrix in workspace. If you wanna
%work with a series of days, go to run_gain_plots and later to
%compile_gain_plots|RubenSep2020


%this code takes the record matrix in day_out.mat and extracts parameter
%space and measurementes (torque, bfields,etc) in day_gain.mat inside

%------------------------------------------------------
%MODIFIED to include Bfiled applied
%Modified to include ALL gauss coeficcients May26,20 Ruben
%Modified to generate GAUSS coeff and not the Blm

%Modified to include many days as well as sinlge day measurements

%DAYS FORMAT EXPLAINED:
% 110212_1 the '_1' means a particular run of the day 110212. For instance 
%you did first a mag ramp, then a reynolds ramp. so 110212_1 is the mag
%ramp and 110212_2 is the rey ramp

%the *out.mat extention is a workspace that includes the raw data aka the
%record matrix from artur's gain_3m_data_chunk code together with t1 t1 tb1
%and tb2. I normally add a variable INFO with description of the run.


%Modified to do gcoeff(31thermaldebiassedprobes - coil signals) dic2020
%previously gains_plot did gcoeff(31thermaldebiassedprobes) - gcoeff(coil signals)

%This is the new gains_plot2 updated from the newest gain_plot feb2022


% function []=gains_plot(run)

%MANY DAYS ------------------------------------------
% masterdir=dir([pwd '/*out.mat']);

% %SINGLE DAY---------------------------------------
% masterdir=dir([pwd '/081022_3out.mat']); 

masterdir=dir([pwd '/' run 'out.mat']); %for function or manual



clearvars i j k
for j=1:length(masterdir)
    folder=masterdir(j).name(1:end-9);
    
    load(masterdir(j).name,'-regexp', '^(?!i|j|masterdir)\w'); %loads everything in masterdirj except everything after ! separated by |
    
    tim=zeros(1,size(record,1));
    ro1=zeros(1,size(record,1));
    re=ones(1,size(record,1))*((2*pi*(1.05^2))/0.7e-6);
    rm=ones(1,size(record,1))*((2*pi*(1.05^2))/0.079);
    out_freq=zeros(1,size(record,1));
    inn_freq=zeros(1,size(record,1));
    torinn=zeros(1,size(record,1));
    powinn=zeros(1,size(record,1));
    powtot=zeros(1,size(record,1));
    g=ones(1,size(record,1))./(((0.7e-6)^2)*(0.45)*927);
%     gginf=zeros(1,size(record,1)); %doing this now in norm_2_gain Feb2022
%     if str2double(folder(5:6)) >= 21
%     fun_G_Re_3M = @(l) (l^1.845416072636967)*0.041601147102921+1.313735481971338e+11;
%     else
%     fun_G_Re_3M = @(l) (l^1.885)*0.003235;
%     end
    % torqout=zeros(1,size(record,1)); no torque out in record yet
    bext=zeros(1,size(record,1));
    
    vrad=zeros(1,size(record,1));
    vphi=zeros(1,size(record,1));
    vrad2=zeros(1,size(record,1));
    vphi2=zeros(1,size(record,1));
    varray=zeros(size(record,1),31); %hall probe from 1-31, in the outer shell
    garray=zeros(size(record,1),24);
    %     bdip = zeros(1,size(record,1));
    %     bcua = zeros(1,size(record,1));
    %
    gl1m0 = zeros(1,size(record,1));
    gl1mc1 =zeros(1,size(record,1));
    gl1ms1 =zeros(1,size(record,1));
    gl1m1 =zeros(1,size(record,1));
    gl2m0 =zeros(1,size(record,1));
    gl2mc1 =zeros(1,size(record,1));
    gl2ms1 =zeros(1,size(record,1));
    gl2m1 =zeros(1,size(record,1));
    gl2mc2 =zeros(1,size(record,1));
    gl2ms2 =zeros(1,size(record,1));
    gl2m2 =zeros(1,size(record,1));
    gl3m0 =zeros(1,size(record,1));
    gl3mc1 =zeros(1,size(record,1));
    gl3ms1 =zeros(1,size(record,1));
    gl3m1 =zeros(1,size(record,1));
    gl3mc2 =zeros(1,size(record,1));
    gl3ms2 =zeros(1,size(record,1));
    gl3m2 =zeros(1,size(record,1));
    gl3mc3 =zeros(1,size(record,1));
    gl3ms3 =zeros(1,size(record,1));
    gl3m3 =zeros(1,size(record,1));
    gl4m0 =zeros(1,size(record,1));
    gl4mc1 =zeros(1,size(record,1));
    gl4ms1 =zeros(1,size(record,1));
    gl4m1 =zeros(1,size(record,1));
    gl4mc2 =zeros(1,size(record,1));
    gl4ms2 =zeros(1,size(record,1));
    gl4m2 =zeros(1,size(record,1));
    gl4mc3 =zeros(1,size(record,1));
    gl4ms3 =zeros(1,size(record,1));
    gl4m3 =zeros(1,size(record,1));
    gl4mc4 =zeros(1,size(record,1));
    gl4ms4 =zeros(1,size(record,1));
    gl4m4 =zeros(1,size(record,1));
    
    %Definning var for coils configuration________________________________%
    
    if str2double(folder(5:6)) >= 15 && str2double(folder(5:6)) < 21  %this is for chosing coils configuration
        var=2;
        if strcmp(folder,'032416') || strcmp(folder,'111715')
            var=3;
        else
        end
    else
        var=1;
    end
    
    if str2double(folder(5:6)) >= 21
        var=7;
    end
    
    %____________________________________________________________________%
    
    
    for i=1:size(record,1)
        tim(i)=record{i,3}{1,1}(1,2);
        ro1(i)=(1/record{i,2}(1));
        re(i)=abs(-record{i,2}(2) + record{i,2}(3))*re(i);
        rm(i)=abs(-record{i,2}(2) + record{i,2}(3))*rm(i);
          
        out_freq(i)=(record{i,2}(1,2));
        inn_freq(i)=(record{i,2}(1,3));
        torinn(i)=mean(record{i,4}{1,2}(:,1)).*1130/1801990;%constant found in Doug codes used to normalized torque data. whERE THE heck comes from, IDK update22 we found where it came from
        g(i)=abs(torinn(i)*g(i));
        %Im also doing abs to
        %compensate rotation direc.
%         gginf(i)=g(i)/fun_G_Re_3M(re(i)); %%doing this now in norm_2_gain Feb2022
        %         torqout(i)=mean(record{i,4}{1,2}(:,1));
        powinn(i)=abs(torinn(i)*inn_freq(i))*2*pi/1000; %in kWatts
        powtot(i)=abs(torinn(i)*(inn_freq(i)-out_freq(i)))*2*pi/1000; %in kWatt
        bext(i)=record{i, 3}{1, 1}(1); %in amps
        if isnan(bext(i)); bext(i) = 0;end %in case the magnets are off it sets nan as zero
        vrad(i)=mean(record{i,3}{1,4}(:,32)); %when i normaliced in norm_gain_plot then scaling from volts to tesla doesn't matter
        vphi(i)=mean(record{i,3}{1,4}(:,33));
        if str2double(folder(5:6)) >= 21
        vrad2(i)=mean(record{i,3}{1,4}(:,34)); %when i normaliced in norm_gain_plot then scaling from volts to tesla doesn't matter
        vphi2(i)=mean(record{i,3}{1,4}(:,35));
        end
        
        %Stablishing coils signal to substract before doing gcoeff ______%
        
            if bext(i)<1  %this is in case of 0 external applied field
                coil = zeros(1,35);
            else
                coil=coils_signal_poly(bext(i),var); %this is a 33 or 35 vector, being the last ones brad and bphi
            end 
            
        
        for k=1:31
            varray(i,k)=mean(record{i,3}{1,4}(:,k)); %hall array on top of each other
        end
        
        
   
        garray(i,:)=gcoeff3m((varray(i,:)-coil(1:31)),probepos); %gauss coeff of each step in parameter space on top of each other
       
        
        gl1m0(i)  =(garray(i,1));
        gl1mc1(i) =(garray(i,2));
        gl1ms1(i) =(garray(i,3));
        gl1m1(i)  =(sqrt(garray(i,2).^2 + garray(i,3).^2));
        gl2m0(i)  =(garray(i,4));
        gl2mc1(i) =(garray(i,5));
        gl2ms1(i) =(garray(i,6));
        gl2m1(i) = (sqrt(garray(i,5).^2 + garray(i,6).^2));
        gl2mc2(i) =(garray(i,7));
        gl2ms2(i) =(garray(i,8));
        gl2m2(i) = (sqrt(garray(i,7).^2 + garray(i,8).^2));
        gl3m0(i)  =(garray(i,9));
        gl3mc1(i) =(garray(i,10));
        gl3ms1(i) =(garray(i,11));
        gl3m1(i)  =(sqrt(garray(i,10).^2 + garray(i,11).^2));
        gl3mc2(i) =(garray(i,12));
        gl3ms2(i) =(garray(i,13));
        gl3m2(i) = (sqrt(garray(i,12).^2 + garray(i,13).^2));
        gl3mc3(i) =(garray(i,14));
        gl3ms3(i) =(garray(i,15));
        gl3m3(i) = (sqrt(garray(i,14).^2 + garray(i,15).^2));
        gl4m0(i)  =(garray(i,16));
        gl4mc1(i) =(garray(i,17));
        gl4ms1(i) =(garray(i,18));
        gl4m1(i) = (sqrt(garray(i,17).^2 + garray(i,18).^2));
        gl4mc2(i) =(garray(i,19));
        gl4ms2(i) =(garray(i,20));
        gl4m2(i) = (sqrt(garray(i,19).^2 + garray(i,20).^2));
        gl4mc3(i) =(garray(i,21));
        gl4ms3(i) =(garray(i,22));
        gl4m3(i) = (sqrt(garray(i,21).^2 + garray(i,22).^2));
        gl4mc4(i) =(garray(i,23));
        gl4ms4(i) =(garray(i,24));
        gl4m4(i) = (sqrt(garray(i,23).^2 + garray(i,24).^2));
     
        
    end
    
    fprintf('Going through file %1.0f of %1.0f.\n',j,numel(masterdir));
%     save([masterdir(j).name(1:end-8) 'amp.mat']);
%     save([masterdir(j).name(1:end-8) '_amp.mat'], '-regexp', '^(?!record|masterdir)\w')
save([masterdir(j).name(1:end-7), 'amp.mat'],'powtot','powinn','out_freq','inn_freq','vphi2','vrad2','vphi','vrad','varray','folder','tim','ro1','re','rm','torinn', 'g','bext','gl1m0','gl1mc1','gl1ms1','gl2m0','gl2mc1','gl2ms1','gl2mc2','gl2ms2','gl3m0','gl3mc1','gl3ms1','gl3mc2','gl3ms2','gl3mc3','gl3ms3','gl4m0','gl4mc1','gl4ms1','gl4mc2','gl4ms2','gl4mc3','gl4ms3','gl4mc4','gl4ms4','gl1m1','gl2m1','gl2m2','gl3m1','gl3m2','gl3m3','gl4m1','gl4m2','gl4m3','gl4m4');
% %
end

clearvars i j k

% end
% '-v7.3'



