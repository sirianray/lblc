% debugging tool for bounce-back-link
% It calculates number of links a lattice point [i1, j1, k1] has associated to drop center [x0, y0, z0]

x0=23; y0=26;
%x0=20; y0=20;
z0=20.5;
r=10;

%i1=12; j1=24; k1=19;
i1=33; j1=31; k1=23;

Lx=40;
Ly=40;
Lz=41;

e=[0 1 -1 0 0 0 0 1 1 1 1 -1 -1 -1 -1];
e=[e; 0 0 0 1 -1 0 0 1 1 -1 -1 1 1 -1 -1];
e=[e; 0 0 0 0 0 1 -1 1 -1 1 -1 1 -1 1 -1];
e=e';

dx=i1-x0;
dy=j1-y0;
dz=k1-z0;
if dx<-Lx/2;
    dx=dx+Lx;
elseif dx>Lx/2
    dx=dx-Lx;
end
if dy<-Ly/2;
    dy=dy+Ly;
elseif dy>Ly/2
    dy=dy-Ly;
end
dr=sqrt(dx^2+dy^2+dz^2);
if dr<=r
    p1=1;
else
    p1=-1;
end

ilink=0;
for i=2:15
    i2=i1+e(i,1);
    j2=j1+e(i,2);
    k2=k1+e(i,3);
    dx=i2-x0;
    dy=j2-y0;
    dz=k2-z0;
    if dx<-Lx/2;
        dx=dx+Lx;
    elseif dx>Lx/2
        dx=dx-Lx;
    end
    if dy<-Ly/2;
        dy=dy+Ly;
    elseif dy>Ly/2
        dy=dy-Ly;
    end
    dr=sqrt(dx^2+dy^2+dz^2);
    if dr<=r
        p2=1;
    else
        p2=-1;
    end
    if p1*p2<0
        ilink=ilink+p1;
    end
end
ilink
