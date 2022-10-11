% THIS IS norm2_gain_plot FOR WHEN USING gains_plot2  ruben FEB22
%__________________________________________________________________
%normalicing brad bphi bcua, bdip using coil_signal code and already
%obtained gain_plots. you
%modified to include ALL gauss coefficients 17April2020
%Modified to create B's from Gauss coeffs coming from run_gain_plots 081120
%also doing norm square of gcoeff

%NOTE: Im not using bias files from record matrix. Instead Im using coils
%signal file from artur as well. I think there are not many discrepancy but
%it's important to remember this. Ruben 081020

% function norm2_gain_plots(run)

%MANY DAYS ------------------------------------------
% masterdir=dir([pwd '/*amp.mat']);

%SINGLE DAY-----------------------------------------
% masterdir=dir([pwd '/081022_3amp.mat']);
masterdir=dir([pwd '/' run 'amp.mat']);

%TYPE OF NORMALIZATION -----------------------------

norma = 0; %If norm 1 then we normalize, if norm 0 then we plot in gauss


clearvars i j
for j=1:length(masterdir) %all %this are the days
    % for k=16
    load(masterdir(j).name, '-regexp', '^(?!record|i|j|masterdir|norma)\w') %LOADING EVERYTHING EXCEPT THE FKING HEAVY RECORD MATRIX
    
    
    gginf=zeros(1,length(bext));
    
    bphi=zeros(1,length(bext));
    brad=zeros(1,length(bext));
    bphi2=zeros(1,length(bext));
    brad2=zeros(1,length(bext));
    barray=zeros(length(bext),31);
    
    bl1m0 = zeros(1,length(bext));
    bl1mc1 =zeros(1,length(bext));
    bl1ms1 =zeros(1,length(bext));
    bl1m1 = zeros(1,length(bext));
    bl2m0 = zeros(1,length(bext));
    bl2mc1 =zeros(1,length(bext));
    bl2ms1 =zeros(1,length(bext));
    bl2m1 = zeros(1,length(bext));
    bl2mc2 =zeros(1,length(bext));
    bl2ms2 =zeros(1,length(bext));
    bl2m2 = zeros(1,length(bext));
    bl3m0 = zeros(1,length(bext));
    bl3mc1 =zeros(1,length(bext));
    bl3ms1 =zeros(1,length(bext));
    bl3m1 = zeros(1,length(bext));
    bl3mc2 =zeros(1,length(bext));
    bl3ms2 =zeros(1,length(bext));
    bl3m2 = zeros(1,length(bext));
    bl3mc3 =zeros(1,length(bext));
    bl3ms3 =zeros(1,length(bext));
    bl3m3 = zeros(1,length(bext));
    bl4m0 =zeros(1,length(bext));
    bl4mc1 =zeros(1,length(bext));
    bl4ms1 =zeros(1,length(bext));
    bl4m1 = zeros(1,length(bext));
    bl4mc2 =zeros(1,length(bext));
    bl4ms2 =zeros(1,length(bext));
    bl4m2 = zeros(1,length(bext));
    bl4mc3 =zeros(1,length(bext));
    bl4ms3 =zeros(1,length(bext));
    bl4m3 = zeros(1,length(bext));
    bl4mc4 =zeros(1,length(bext));
    bl4ms4 =zeros(1,length(bext));
    bl4m4 = zeros(1,length(bext));
    
    if str2double(folder(5:6)) >= 15 && str2double(folder(5:6)) < 21  %this is for chosing coils configuration
        var=2;
        if strcmp(folder,'032416') || strcmp(folder,'111715') %selecting cuadrupoes modes
            var=3;
        else
        end
    elseif str2double(folder(5:6)) <= 14
        var=1;
    elseif str2double(folder(5:6)) >= 21
        var=7;
    end
    
    for i=1:length(bext) %this is every parameter of every particular run/day
        
        
        if norma == 1
            
            
            if bext(i)<1  %this is in case of 0 external applied field
                coil=coils_signal_poly(20,var);
                fprintf('Element path = %1.0f normalized with respect to 20 amps in the coils.\n',i)
                
            else
                coil=coils_signal_poly(bext(i),var); %this is a 33 vector, being the last ones brad and bphi
            end
            %                 coil_coeff=gcoeff3m(coil(1:31),probepos);%this is a 24 vector being the 24 gauss coeff l=4
            coil_coeff = zeros(1,24); %doing this if gainplot2 was used
            
            norm_fac = coil(32);
        end
        
        if norma == 0
            %           fprintf('Bfield in Gauss.\n')
            norm_fac = 0.0315;
            if bext(i)<1
                coil_coeff = zeros(1,24); %this was first attemp of plotting without substracting coils
                coil = zeros(1,35);
            else
                
                coil=coils_signal_poly(bext(i),var); %this is a 33 vector, being the last ones brad and bphi
                coil_coeff = zeros(1,24); %doing this if gainplot2 was used
                %             coil_coeff=gcoeff3m(coil(1:31),probepos);%this is a 24 vector being the 24 gauss coeff l=4
                
            end
            
        else
        end
        
        
        if str2double(folder(5:6)) >= 21
            fun_G_Re_3M = @(l) (l^1.845416072636967)*0.041601147102921+1.313735481971338e+11;
        else
%             fun_G_Re_3M = @(l) (l^1.885)*0.003235;  % This is the old one used for smooth before Abril 2022. I was not taken into acount friction
                fun_G_Re_3M = @(l) (l^1.83)*0.009 + 2.366e10;
        end
        
        gginf(i)=g(i)/fun_G_Re_3M(re(i));
        
        
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
        %THIS IS THE IMPORTANT FIN________________________________________________________________%
        
        brad(i)=(vrad(i)-coil(32))/norm_fac;
        bphi(i)=(vphi(i)-coil(33))/norm_fac_phi;
        
        if str2double(folder(5:6)) >= 21
            brad2(i)=(vrad2(i)-coil(34))/norm_fac;
            bphi2(i)=(vphi2(i)-coil(35))/norm_fac;
        end
        
        
        for k=1:31
            barray(i,k)=(varray(i,k)-coil(k))/norm_fac; %hall array on top of each other
        end
        
        bl1m0(i)  =(2*0.95^3*(gl1m0(i)-coil_coeff(1))/norm_fac);
        bl1mc1(i) =2*0.95^3*(gl1mc1(i)-coil_coeff(2))/norm_fac;
        bl1ms1(i) =2*0.95^3*(gl1ms1(i)-coil_coeff(3))/norm_fac;
        bl1m1(i)= 2*0.95^3*(gl1m1(i)-sqrt(coil_coeff(2)^2+coil_coeff(3)^2))/norm_fac; %this sqrt(coil_coeff...) its overdoing it. Average or Higher would be enough.
        bl2m0(i)  =(6*0.95^4*(gl2m0(i)-coil_coeff(4))/norm_fac);
        bl2mc1(i) =6*0.95^4*(gl2mc1(i)-coil_coeff(5))/norm_fac;
        bl2ms1(i) =6*0.95^4*(gl2ms1(i)-coil_coeff(6))/norm_fac;
        bl2m1(i)=6*0.95^4*(gl2m1(i)-sqrt(coil_coeff(5)^2+coil_coeff(6)^2))/norm_fac;
        bl2mc2(i) =6*0.95^4*(gl2mc2(i)-coil_coeff(7))/norm_fac;
        bl2ms2(i) =6*0.95^4*(gl2ms2(i)-coil_coeff(8))/norm_fac;
        bl2m2(i)=6*0.95^4*(gl2m2(i)-sqrt(coil_coeff(7)^2+coil_coeff(8)^2))/norm_fac;
        bl3m0(i)  =(12*0.95^5*(gl3m0(i)-coil_coeff(9))/norm_fac);
        bl3mc1(i) =12*0.95^5*(gl3mc1(i)-coil_coeff(10))/norm_fac;
        bl3ms1(i) =12*0.95^5*(gl3ms1(i)-coil_coeff(11))/norm_fac;
        bl3m1(i)=12*0.95^5*(gl3m1(i)-sqrt(coil_coeff(10)^2+coil_coeff(11)^2))/norm_fac;
        bl3mc2(i) =12*0.95^5*(gl3mc2(i)-coil_coeff(12))/norm_fac;
        bl3ms2(i) =12*0.95^5*(gl3ms2(i)-coil_coeff(13))/norm_fac;
        bl3m2(i)=12*0.95^5*(gl3m2(i)-sqrt(coil_coeff(13)^2+coil_coeff(12)^2))/norm_fac;
        bl3mc3(i) =12*0.95^5*(gl3mc3(i)-coil_coeff(14))/norm_fac;
        bl3ms3(i) =12*0.95^5*(gl3ms3(i)-coil_coeff(15))/norm_fac;
        bl3m3(i)=12*0.95^5*(gl3m3(i)-sqrt(coil_coeff(14)^2+coil_coeff(15)^2))/norm_fac;
        bl4m0(i)  =(20*0.95^6*(gl4m0(i)-coil_coeff(16))/norm_fac);
        bl4mc1(i) =20*0.95^6*(gl4mc1(i)-coil_coeff(17))/norm_fac;
        bl4ms1(i) =20*0.95^6*(gl4ms1(i)-coil_coeff(18))/norm_fac;
        bl4m1(i) = 20*0.95^6*(gl4m1(i)-sqrt(coil_coeff(17)^2+coil_coeff(18)^2))/norm_fac;
        bl4mc2(i) =20*0.95^6*(gl4mc2(i)-coil_coeff(19))/norm_fac;
        bl4ms2(i) =20*0.95^6*(gl4ms2(i)-coil_coeff(20))/norm_fac;
        bl4m2(i)=20*0.95^6*(gl4m2(i)-sqrt(coil_coeff(19)^2+coil_coeff(20)^2))/norm_fac;
        bl4mc3(i) =20*0.95^6*(gl4mc3(i)-coil_coeff(21))/norm_fac;
        bl4ms3(i) =20*0.95^6*(gl4ms3(i)-coil_coeff(22))/norm_fac;
        bl4m3(i)=20*0.95^6*(gl4m3(i)-sqrt(coil_coeff(21)^2+coil_coeff(22)^2))/norm_fac;
        bl4mc4(i) =20*0.95^6*(gl4mc4(i)-coil_coeff(23))/norm_fac;
        bl4ms4(i) =20*0.95^6*(gl4ms4(i)-coil_coeff(24))/norm_fac;
        bl4m4(i)=20*0.95^6*(gl4m4(i)-sqrt(coil_coeff(23)^2+coil_coeff(24)^2))/norm_fac;
    end
    
    
    %  save(masterdir(j).name, '-regexp', '^(?!record|i|j|masterdir)\w');
    % %     save(masterdir(j).name(1:end-8), 'folder','tim','ro1','re','rm','torinn','gginf', 'g','bext','brad','bphi','bdip','bcua','brad','bphi','bdip','bcua','bl1m0','bl1mc1','bl1ms1','bl2m0','bl2mc1','bl2ms1','bl2mc2','bl2ms2','bl3m0','bl3mc1','bl3ms1','bl3mc2','bl3ms2','bl3mc3','bl3ms3','bl4m0','bl4mc1','bl4ms1','bl4mc2','bl4ms2','bl4mc3','bl4ms3','bl4mc4','bl4ms4','bl1m1','bl2m1','bl2m2','bl3m1','bl3m2','bl3m3','bl4m1','bl4m2','bl4m3','bl4m4');
    %         clear tim ro1 re rm torinn gginf g brad bphi bdip bcua  bl1m0 bl1mc1 bl1ms1 bl2m0 bl2mc1 bl1ms1 bl2mc2 bl2ms2 bl3m0 bl3mc1 bl3ms1 bl3mc2 bl3ms2 bl3mc3 bl3ms3 bl4m0 bl4mc1 bl4ms1 bl4mc2 bl4ms2 bl4mc3 bl4ms3 bl4mc4 bl4ms4 bl1m1 bl2m1 bl2m2 bl3m1 blm3m2 bl3m3 bl4m1 bl4m2 bl4m3 bl4m4
    save([masterdir(j).name(1:end-7), 'amp.mat'],'powtot','powinn','out_freq','inn_freq','norma','coil','bphi2','brad2','vphi2','vrad2','vphi','vrad','var','varray','folder','tim','ro1','re','rm','torinn','gginf', 'g','bext','brad','bphi','barray','bl1m0','bl1mc1','bl1ms1','bl2m0','bl2mc1','bl2ms1','bl2mc2','bl2ms2','bl3m0','bl3mc1','bl3ms1','bl3mc2','bl3ms2','bl3mc3','bl3ms3','bl4m0','bl4mc1','bl4ms1','bl4mc2','bl4ms2','bl4mc3','bl4ms3','bl4mc4','bl4ms4','bl1m1','bl2m1','bl2m2','bl3m1','bl3m2','bl3m3','bl4m1','bl4m2','bl4m3','bl4m4','gl1m0','gl1mc1','gl1ms1','gl2m0','gl2mc1','gl2ms1','gl2mc2','gl2ms2','gl3m0','gl3mc1','gl3ms1','gl3mc2','gl3ms2','gl3mc3','gl3ms3','gl4m0','gl4mc1','gl4ms1','gl4mc2','gl4ms2','gl4mc3','gl4ms3','gl4mc4','gl4ms4','gl1m1','gl2m1','gl2m2','gl3m1','gl3m2','gl3m3','gl4m1','gl4m2','gl4m3','gl4m4');
    clear k j i
end

% end %function end
