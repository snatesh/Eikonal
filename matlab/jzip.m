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
m = 14; n = m+1;
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
% vandermonde under (a,b,c), (a+1,b,c) etc.
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
% make coeffs of fref under (a,b,c)
cfref_abc = V_abc'*FrefW;
disp(norm(V_abc*cfref_abc-Fref)/normFref);

% reference tri
Rv = [0,1,0];
Sv = [0,0,1];

img = imread("mathias-reding-vU9-VO-4Nk0-unsplash.jpg");
img = rgb2gray(img);
imgsz = size(img);
nXpix = imgsz(1); nYpix = imgsz(2);
[Xpix,Ypix] = meshgrid(1:nXpix,1:nYpix); 

nSamp = 10;
[X,Y] = meshgrid(linspace(1,nXpix,nSamp), linspace(1,nYpix,nSamp));
X = X(:); Y = Y(:);
DT = delaunayTriangulation(X,Y);
triInd = DT.pointLocation(Xpix(:),Ypix(:));
T = DT.ConnectivityList;
figure(1);
triplot(T,X,Y); hold on;
imagesc(img','AlphaData',0.5)
nTri = length(T);
interpolator = @(x,y) interp2(Xpix,Ypix,double(img)',x,y,'makima');
Xpix = Xpix(:); Ypix = Ypix(:);
Imgapprox = zeros(size(img'));
for j = 1:nTri

    Xe = [X(T(j,:))';Y(T(j,:))'];
    Ixe = IncidenceMatrix(Xe);
    XYe = (Ixe * [R,S]' + Xe(:,1))';
    imginterp = interpolator(XYe(:,1),XYe(:,2));
    cimg = V_abc' * (imginterp .* W);
    % find pixels inside current triangle
    xpix = Xpix(triInd==j);
    ypix = Ypix(triInd==j);
    % map those pixels to reference triangle
    XYpix = Ixe \ [xpix'-Xe(1,1);ypix'-Xe(2,1)];
    % evaluate vandermonde at pixels in reference
    Vabcpix = jPoly_tri(XYpix(1,:)',XYpix(2,:)',H_abc,n-1,a,b,c);
    % evaluate approximation to image within triangle
    imgapprox = Vabcpix*cimg;
    Imgapprox(triInd==j) = round(imgapprox);
    % get actual pixel values within triangle
    imgactual = interpolator(xpix,ypix);
    disp(norm(round(imgapprox)-imgactual)/norm(imgactual));
    %plot(XYe(:,1),XYe(:,2),'k.');
    %plot(XYpix(1,:)',XYpix(2,:)','g.');
    %drawnow; pause(0.1);
    disp(j)
end


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
