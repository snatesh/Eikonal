clear all; close all; clc;
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
  default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
  set(groot, default_name,'latex');
end
set(groot, 'defaultLegendFontSize',30)
set(groot, 'defaultAxesFontSize',30)
set(groot, 'defaultLineLineWidth',3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% jacobi poly params
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
load './triquadLeg_17_28.mat';
R = Zk(1:N); S = Zk(N+1:2*N); W = Zk(2*N+1:3*N);
% weight function for (a+1,b+1,c+1)
w_a1b1c1 = gamma(a+1+b+1+c+1+3/2)/(gamma(a+1+1/2)*gamma(b+1+1/2)*gamma(c+1+1/2));
wa1b1c1 = @(x,y) x.^(a+1-1/2).*y.^(b+1-1/2).*(1-x-y).^(c+1-1/2)*w_a1b1c1/2;
% test function
f = @(x,y) (x+y).^10; 
% eval test function on quadrature nodes
Fref = f(R,S); FrefW = Fref.*W; normFref = norm(Fref);
m = 17; n = m+1;
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
% vandermonde under (a,b,c), (a+1,b,c) etc.
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);

% make coeffs of fref under (a,b,c), (a+1,b,c) etc.
cfref_abc = V_abc'*FrefW;
disp(norm(V_abc*cfref_abc-Fref)/normFref);

% reference tri
Rv = [0,1,0];
Sv = [0,0,1];
Xeref = [Rv;Sv];

% general tri
Xv = [2.,4.,3.];
Yv = [1.,1.,2.];
Xe = [Xv;Yv];

Ixe = IncidenceMatrix(Xe);
XY = Ixe * [R,S]' + Xe(:,1);
XY = XY';
X = XY(:,1); Y = XY(:,2);
F = f(X,Y); FW = F.*W; normF = norm(F);
% make coeffs of fref under (a,b,c), (a+1,b,c) etc.
cf_abc = V_abc'*FW;
disp(norm(V_abc*cf_abc-F)/normF);

%scatter3(R,S,F); hold on;
%scatter3(X,Y,F);
img = imread("mathias-reding-vU9-VO-4Nk0-unsplash.jpg");
img = rgb2gray(img);
imgsz = size(img);
Xpix = imgsz(1); Ypix = imgsz(2);
Xp = 1:Xpix; Yp = 1:Ypix;
[XX,YY] = meshgrid(Yp,Xp);
Xv1 = [Xp(1),Xp(end),Xp(1)];
Yv1 = [Yp(1),Yp(1),Yp(end)];
Xv2 = [Xp(1),Xp(end),Xp(end)];
Yv2 = [Yp(end),Yp(1),Yp(end)];
Xe1 = [Xv1;Yv1];
Xe2 = [Xv2;Yv2];
plot_tri(Xv1,Yv1,'k-');
plot_tri(Xv2,Yv2,'b-');

Ixe1 = IncidenceMatrix(Xe1);
Ixe2 = IncidenceMatrix(Xe2);
XY1 = (Ixe1 * [R,S]' + Xe1(:,1))';
XY2 = (Ixe2 * [R,S]' + Xe2(:,1))';
plot(XY1(:,1),XY1(:,2),'k.');
plot(XY2(:,1),XY2(:,2),'b.');
X = [XY1(:,1);XY2(:,1)];
Y = [XY1(:,2);XY2(:,2)];
plot(X,Y,'go')
imginterp = interp2(Xp,Yp,double(img)',X,Y,'makima');
cimg1 = V_abc'*(imginterp(1:N).*W);
cimg2 = V_abc'*(imginterp(N+1:end).*W);
hold off;
scatter3(X,Y,imginterp,'rp'); hold on;
scatter3(XY1(:,1),XY1(:,2),V_abc*cimg1,'ko')
scatter3(XY2(:,1),XY2(:,2),V_abc*cimg2,'bo')

%%

img = imread("mathias-reding-vU9-VO-4Nk0-unsplash.jpg");
imgflat = double(img(:)); 
N = length(imgflat);

dctimg = dct(imgflat);
dctimg_sort = dct(sort(imgflat));
semilogy(abs(dctimg),'bo-'); hold on;
semilogy(abs(dctimg_sort),'rs-');


%%
plot_tri(Rv,Sv,'k-');
plot_tri(Xv,Yv,'b-');


plot(R,S,'k.');
plot(X,Y,'b.');

function Ixe = IncidenceMatrix(Xe)

Ixe = [Xe(:,2)-Xe(:,1), Xe(:,3)-Xe(:,1)];

end

function plot_tri(Xv,Yv,mrk)

e1x = [Xv(1),Xv(2)];
e1y = [Yv(1),Yv(2)];
e2x = [Xv(2),Xv(3)];
e2y = [Yv(2),Yv(3)];
e3x = [Xv(3),Xv(1)];
e3y = [Yv(3),Yv(1)];

plot(e1x,e1y,mrk); hold on;
plot(e2x,e2y,mrk);
plot(e3x,e3y,mrk);

end
