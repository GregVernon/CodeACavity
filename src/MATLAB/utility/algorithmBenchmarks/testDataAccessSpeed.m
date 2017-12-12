clear
N = 100;
x = linspace(0,1,N);
y = linspace(0,1,N);
[xx,yy] = meshgrid(x,y);

ii = 1:N;
jj = 1:N;
[ii,jj] = meshgrid(ii,jj);

DATA.xx = xx;
DATA.yy = yy;
DATA.ii = ii;
DATA.jj = jj;

zz = DATA.xx(DATA.jj, DATA.ii) + DATA.yy(DATA.jj,DATA.ii);
ZZ = xx(jj,ii) + yy(jj,ii);

k = zz+2;
K = ZZ+2;
