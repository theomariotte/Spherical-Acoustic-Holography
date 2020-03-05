function dessin_13HP(x0,y0)
load caissonHP13
% Démarche : chercher les coordonnées du HP central (x0,y0)
% dans le repère de la grille scannée
% par défaut : x0=0.33, y0=0.33
% centre de gravité du rectangle
xcentre=0.355;
ycentre=0.355;


if (nargin==0)
    x0=0.33;
    y0=0.33;
else
    x0=x0-xcentre;
    y0=y0-ycentre;
end

nbrec=size(rectangles,1);
hold on
p1=[rectangles(:,1)-rectangles(:,3)/2 rectangles(:,2)+rectangles(:,4)/2]+ones(nbrec,1)*[x0 y0];
p2=[rectangles(:,1)+rectangles(:,3)/2 rectangles(:,2)+rectangles(:,4)/2]+ones(nbrec,1)*[x0 y0];
p3=[rectangles(:,1)+rectangles(:,3)/2 rectangles(:,2)-rectangles(:,4)/2]+ones(nbrec,1)*[x0 y0];
p4=[rectangles(:,1)-rectangles(:,3)/2 rectangles(:,2)-rectangles(:,4)/2]+ones(nbrec,1)*[x0 y0];
for i=1:nbrec
    Pts=[p1(i,:);p2(i,:);p3(i,:);p4(i,:);p1(i,:)];
    plot(Pts(:,1),Pts(:,2),'LineWidth',linewidth,'color',linecolor)
end
% if (strcmp(get(gca,'xdir'),'reverse'))
%     circles(3:end,1)=-circles(3:end,1);
% end
theta=linspace(0,2*pi,128);
ctheta=cos(theta);
stheta=sin(theta);
nbcer=size(circles,1);
for i=1:nbcer
    x=circles(i,3)/2*ctheta+circles(i,1)+x0;
    y=circles(i,3)/2*stheta+circles(i,2)+y0;
    plot(x,y,'LineWidth',linewidth,'color',linecolor)
end
%axis equal
axis square
hold off