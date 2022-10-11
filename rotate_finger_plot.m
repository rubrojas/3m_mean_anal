%Rotation Matrix plOT
figure(39);clf;
run='081722_1';

cd(['/Users/rubenrojas/Documents/3M/Amplification/3m_folder_files/' run(1:6)]);
load(['/Users/rubenrojas/Documents/3M/Amplification/3m_folder_files/' run(1:6) '/' run 'amp.mat']);
update = 1;
x = rm;
bdif = zeros(36,1); %this has the length of theta
bfix = zeros(length(brad),2);
theta = 10:10:360;
thetamin = zeros(length(brad),1);
for i=1:length(brad)
for j=1:36
b2 = [brad(i) bphi(i)];
b = [brad2(i) bphi2(i)];
R = [cosd(theta(j)) -sind(theta(j)); sind(theta(j)) cosd(theta(j))];
bnew = b2*R;
bdif(j)=norm(b-bnew);
% scatter(theta(j),bdif(j),'filled','k');

end
[m,ind] = min(bdif);
thetamin(i) = theta(ind);
% thetamin(i) = 50;
bfix(i,:) = b2*[cosd(thetamin(i)) -sind(thetamin(i)); sind(thetamin(i)) cosd(thetamin(i))];
end

scatter(1:length(brad),thetamin,'filled','k');


figure(40);clf;subplot(1,2,1);plot(x,bfix(:,1),'b');hold on;plot(x,brad,'b');plot(x,bfix(:,2),'r');plot(x,bphi,'r');
subplot(1,2,2);plot(x,brad2,'b');hold on;plot(x,brad,'b');plot(x,bphi,'r');plot(x,bphi2,'r');

 brad2o=brad2;
 bphi2o=bphi2;
 brad2=bfix(:,1);
 bphi2=bfix(:,2);
 

 if update == 1
      fprintf('Press any key to update brad2 and bphi2...\n'); 
     pause;
     save([run 'amp.mat'],'brad2','bphi2','brad2o','bphi2o','-append');
 end
 
%  %undo
%  brad2=brad2o;
% bphi2=bphi2o;
 

