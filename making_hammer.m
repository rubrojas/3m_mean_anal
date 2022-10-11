% Perevalov Sep 2019
% creates a picture with Aitoff projection using known gauss coefficients
% REQUIRES m_map package https://www.eoas.ubc.ca/~rich/map.html#download

%Ruben 2021 making this a funcntion
% adding lim to fix the color scale

function [] = making_hammer(fake_gauss,lim)
%% Here are the gauss coeeficients array for test
% fake_gauss = zeros(1,24);
% fake_gauss(1) = 3.3;   % this one just adds something non-zero
% fake_gauss(9) = -4;
% fake_gauss(16) = -2.2;
% coils=coils_signal(200,3);
% fake_gauss=gcoeff3m(coils(1:31),probepos);
% fake_gauss = bcoeff(600303,:);
%% This is for the actual data
% construct g_array from data vectors on workspace

% m=0      m=1(cos),m=-1(sin)      m=2(cos),m=-2(sin)      m=3(cos),m=-3(sin)      m=4(cos),m=-4(sin)             
% gcarray_cs= {bl1m0,{bl1mc1,bl1ms1},       {NaN,NaN},              {NaN,NaN},              {NaN,NaN};        %l=1
%                     bl2m0,{bl2mc1,bl2ms1},{bl2mc2,bl2ms2},       {NaN,NaN},              {NaN,NaN};        %l=2
%                     bl3m0,{bl3mc1,bl3ms1},{bl3mc2,bl3ms2},{bl3mc3,bl3ms3},       {NaN,NaN};        %l=3
%                     bl4m0,{bl4mc1,bl4ms1},{bl4mc2,bl4ms2},{bl4mc3,bl4ms3},{bl4mc4,bl4ms4}};%l=4


% fake_gauss=horzcat(bl1m0'./2,bl1mc1'./2,bl1ms1'./2,bl2m0'./6,bl2mc1'./6,bl2ms1'./6,bl2mc2'./6,bl2ms2'./6,bl3m0'./12,bl3mc1'./12,bl3ms1'./12,bl3mc2'./12,bl3ms2'./12,bl3mc3'./12,bl3ms3'./12,bl4m0'./20,bl4mc1'./20,bl4ms1'./20,bl4mc2'./20,bl4ms2'./20,bl4mc3'./20,bl4ms3'./20,bl4mc4'./20,bl4ms4'./20);
% fake_gauss=coils_signal(100,1);
%  fake_gauss = record{1,3}{1,6}(10000,:);
 
 
 
% creating coordinates to plot onto
phi = (0:0.05:2*pi);
theta=(0:0.05:pi);
phi_deg = phi.*180/pi-180;
theta_deg = theta.*180/pi-90;

grid_coords = zeros(length(phi)*length(theta),3);

% making coordinates of the imaginary probes at these locations
for phi_i = 1:length(phi)
    for theta_i = 1:length(theta)
        k = length(theta)*(phi_i-1)+theta_i;
        grid_coords(k,1) = 1;
        grid_coords(k,2) = theta(theta_i); % theta
        grid_coords(k,3) = phi(phi_i);
    end
end

% returns the values on these imaginary probes based on the g's
% B_array = gauss2_hall(fake_gauss(19,:),grid_coords);
B_array = gauss2_hall(fake_gauss,grid_coords);




% converting from 1D array into 2D 
B_map = zeros(length(phi),length(theta));
for phi_i = 1:length(phi)
    for theta_i = 1:length(theta)
        k = length(theta)*(phi_i-1)+theta_i;
        B_map(phi_i,theta_i) = B_array(k);        
    end
end




% Plotting
DanMap = [49 26 0; 73 39 0;102 54 0;146 78 0 ;174 94 0;202 110 0 ;222 136 0;237 171 0 ;251 222 0 ;255 255 0 ;0 255 255;0 200 232 ;0 149 202;0 118 182;0 98 167;0 68 141;0 50 102;0 23 75;0 0 51;0 6 27];DanMap = (DanMap./255);
fz=30;    
hammer_fig = figure(6);
if ~exist('lim','var')
lim=max(max(B_map));
else
end
lim=1.5;
m_proj('Hammer-Aitoff','lat',[-90 90], 'lon',[-180 180])
m_contourf(phi_deg,theta_deg,B_map','EdgeColor','None')
caxis([-lim lim]);
% caxis([-3 3]);
c=colorbar;
yl=ylabel(c,'$B ({\rm G})$','FontName','Times New Roman','Interpreter','latex','FontSize',fz);
set(yl, 'position',[2.5 0 0]);

%  colormap(DanMap);
colormap(flipud(DanMap));
% colormap(flipud(brewermap(1000,'RdBu')));
box on
set(gca,'linewidth',1,'FontName','Times New Roman','FontSize',fz);
set(gca,'linewidth',1,'Yticklabel',[],'Xticklabel',[]);
set(gcf,'color','w');
pbaspect([2 1 1]);
set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});






end

