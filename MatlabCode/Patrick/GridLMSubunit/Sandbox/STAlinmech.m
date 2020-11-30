% Testing predictions from singular mechanism 

global GLMP
sub = 2;

Lcc = GLMP.subunit{sub}.Lcc;
Mcc = GLMP.subunit{sub}.Mcc;
nsp = GLMP.subunit{sub}.nspikes;

X = [Lcc Mcc];
M = sqrtm(inv(X'*X));

figure(1); clf; hold on;
plot3(Lcc,Mcc,nsp,'K*')
xlabel('Lcc'); ylabel('Mcc'); zlabel('Nsp')

XX = X * M;

figure(2); clf; hold on;
plot(XX(:,1),XX(:,2),'k*')
axis equal tight
box on


