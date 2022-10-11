%RUBEN 2021 Updated from Artur Maiking_aitoff
%Now I want to also make a movie of the evulution of the Bfield at a
%particular parameter space run to see the structure

% creates a picture with Aitoff projection using known gauss coefficients
% REQUIRES m_map package https://www.eoas.ubc.ca/~rich/map.html#download

%Preeliminares for run day------------------------------------------------------
%MAKE SURE YOU LOAD THE RECORD MATRIX YOU WANT

% run = '060414_1'; %run day

%this is for selecting the point within each run
path = 25;

%Name of the gif figure
name_im = [run '_' num2str(path) '.gif'];
% name_im = 'another.gif';

%chooses the porcentaje of the time section to be taken into account. 1/2
%means the upper half of the date etc, to remove transient
%SET TIME RATIO to 0 to manually choose time_i (starting time of the movie
%in ssm's)
time_ratio = 1;
% time_i = 40982;

%Choose 1 if wish to ALSO GENERATE MOVIE or zero if not --------------------------------------------
MOVIE = 0;

%Define interval of movie in seconds
sec = 0.1;
inte = round(256*sec);

%define Duration of the movie in frames
dura = 60;

%TYPE OF NORMALIZATION -----------------------------

norma = 0; %If norm 1 then we normalize, if norm 0 then we plot in gauss

sz=70; %scatter plots size
fz=24; %font size
lz=14; %legend size
pz=2; %plots size


ampdir = dir([pwd '/' run  'amp.mat']); %seems a directory but its only 1 run
outdir = dir([pwd '/' run  'out.mat']);

clearvars('rm','re','ro1','bext','var','folder','INFO','tim')
load(ampdir.name,'rm','re','ro1','bext','var','folder','INFO','tim');

if exist('INFO','var')
    fprintf([INFO '\n'])
end



if ~exist('record','var')  %checking if record matrix which is heavy AF is already loaded in workspace
    fprintf('Loading record cell...\n')
    load(outdir.name,'record');
    
elseif ~strcmp(record{1,1},run(1:end-2))
    fprintf('Loading record cell...\n')
    load(outdir.name,'record');
    % else
    %     fprintf('Using existing record matrix.\n')
    
end




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

% definning color map
cmap=0.8.*[1,0,0.852000000000000;1,0,0;1,0.396000000000000,0;1,0.792000000000000,0;0.812000000000000,1,0;0.416000000000000,1,0;0,1,0.772000000000000;0,0.436000000000000,1;0,0.0400000000000000,1;0.356000000000000,0,1;0.752000000000001,0,1];

%defining torque funtion for 3m
fun_G_Re_3M = @(l) (l^1.885)*0.003235;

clearvars time brad_t bphi_t coil

%Preeliminares for each point in path------------------------------------------------------------------
for k=path
    
    
    
    if str2double(folder(5:6)) >= 15 %this is for chosing coils configuration
        var=2;
        if strcmp(folder,'032416') || strcmp(folder,'111715')
            var=3;
        else
        end
    else
        var=1;
    end
    
    
    
    if norma == 1
        if  str2double(folder(5:6)) >= 21
            norm_fac = coil(34);
        else
            norm_fac = coil(32);
        end
        
        if bext(k)<1  %this is in case of 0 external applied field
            coil=coils_signal_poly(20,var);
            fprintf('Plotting element path = %1.0f normalized with respect to 20 amps in the coils.\n',k)
        else
            coil=coils_signal_poly(bext(k),var); %this is a 33 vector, being the last ones brad and bphi
        end
        coil_coeff=gcoeff3m(coil(1:31),probepos);%this is a 24 vector being the 24 gauss coeff l=4
        
        
    elseif norma == 0
        fprintf('Plotting Bfield in Gauss.\n')
        norm_fac = 0.0315;
        if bext(k)<1
            coil_coeff = zeros(1,24); %this was first attemp of plotting without substracting coils
            coil = zeros(1,34);
        else
        
        coil=coils_signal_poly(bext(k),var); %this is a 33 vector, being the last ones brad and bphi
        coil_coeff=gcoeff3m(coil(1:31),probepos);%this is a 24 vector being the 24 gauss coeff l=4
        end
        
    else
    end
    
    
    %CHANGE time_ratio accordingly to get rid of transient time. Set to 1
    %to take the whole section of data for each parameter of the run.Set to
    %1/2 to choose second half of the parameter ru. etc.
    
    time = record{k,3}{1,2};
    time_k = length(record{k,3}{1,2});
    if time_ratio ~= 0
        time_i = time_k-round(time_k*time_ratio)+1;
    end
    time = time(time_i:end);
    %     time = time-min(time)+1; %setting t1=0
    %     brad_t = record{k,3}{1,4}(time_i:end,32)./coil(32);
    %     bphi_t = record{k,3}{1,4}(time_i:end,33)./coil(32);
    %     brad2_t = record{k,3}{1,4}(time_i:end,34)./coil(32);
    %     bphi2_t = record{k,3}{1,4}(time_i:end,35)./coil(32);
    %     beq1_t = record{k,3}{1,4}(time_i:end,12)./coil(32);%radial probes from array
    %     beq2_t = record{k,3}{1,4}(time_i:end,13)./coil(32);%radial probes from array
    %     beq3_t = record{k,3}{1,4}(time_i:end,14)./coil(32);%radial probes from array
    %
    
    %     if ~exist('gl1m0_t,'var')
    
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
    %     end
    
    
    %Normalizing Time from sec to outer_rot_rate. NOTE:Since we are setting
    %time t1=0 the PHASE IS UNKNOWN
    
    out_freq=abs(record{k,2}(1,2)); %outer frequency in Hz for normalice
    if out_freq <= 0.001
        out_freq = 1;
    end
    inn_freq=abs(record{k,2}(1,3)); %inner frequency for normalice
    
    
    time = time*out_freq;
    
    
    fake_gauss=horzcat(bl1m0_t./2,bl1mc1_t./2,bl1ms1_t./2,bl2m0_t./6,bl2mc1_t./6,bl2ms1_t./6,bl2mc2_t./6,bl2ms2_t./6,bl3m0_t./12,bl3mc1_t./12,bl3ms1_t./12,bl3mc2_t./12,bl3ms2_t./12,bl3mc3_t./12,bl3ms3_t./12,bl4m0_t./20,bl4mc1_t./20,bl4ms1_t./20,bl4mc2_t./20,bl4ms2_t./20,bl4mc3_t./20,bl4ms3_t./20,bl4mc4_t./20,bl4ms4_t./20);
    making_hammer(fake_gauss(1,:))
    axis tight
    set(gca,'nextplot','replacechildren','visible','off')
    f = getframe;
    [im,map] = rgb2ind(f.cdata,256,'nodither');
    im(1,1,1,dura) = 0;
    l = 0;
    for t = 1:inte:(dura*inte)
        l = l + 1;
        making_hammer(fake_gauss(t,:))
        f = getframe;
        im(:,:,1,l) = rgb2ind(f.cdata,map,'nodither');
    end
    if MOVIE
        imwrite(im,map,name_im,'DelayTime',0.5,'LoopCount',inf) %g443800
    end
    %      clearvars
    
    
    
    %Now we make the video
    
    
    
    
    
    
    
    
end
