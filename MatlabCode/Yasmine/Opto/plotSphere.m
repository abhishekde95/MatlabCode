function plotSphere(x,y,z,r,color)

[X,Y,Z] = sphere;

for n = 1:length(x)
    surf(X*r+x(n),Y*r+y(n),Z*r+z(n),'EdgeColor','none','FaceColor',color);
end