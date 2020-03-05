fig=10;
mic=load('data/3Dcam36Officiel.txt','-ascii');
M=36;
r=0.15;

x=mic(:,1);
y=mic(:,2);
z=mic(:,3);

figure(fig)
plot(x,y,'+')
axis square
text(x,y,num2str([1:M]'))
xlabel('x')
ylabel('y')

figure(fig+1)
plot(y,z,'+')
axis square
text(y,z,num2str([1:M]'))
xlabel('y')
ylabel('z')

figure(fig+2)
plot(x,z,'+')
axis square
text(x,z,num2str([1:M]'))
xlabel('x')
ylabel('z')


[xx,yy,zz]=sphere;

figure(fig+3)
mesh(r*xx,r*yy,r*zz)
hold on
plot3(x,y,z,'+r')
xlabel('x')
ylabel('y')
zlabel('z')
axis square
text(x,y,z,num2str([1:M]'))
hold off
axis equal