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
%load './triquadLeg_17_28.mat';
S = importdata("../bin/ytri_N325_n24_M465_m29.txt"); N = length(S);
R = importdata("../bin/xtri_N325_n24_M465_m29.txt");
W = importdata("../bin/wtri_N325_n24_M465_m29.txt");
%R = Zk(1:N); S = Zk(N+1:2*N); W = Zk(2*N+1:3*N);
% test function
f = @(x,y)(x+y).^14
df = @(x,y) 14*(x+y).^13;
% eval test function on quadrature nodes
Fref = f(R,S); FrefW = Fref.*W; normFref = norm(Fref);
m = 14; n = m+1;
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc1 = structure_factors_tri(n+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(n+1,a,b+1,c+1);
% vandermonde under (a,b,c), (a+1,b,c) etc.
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
V_a1bc1 = jPoly_tri(R,S,H_a1bc1,n-1,a+1,b,c+1);
V_ab1c1 = jPoly_tri(R,S,H_ab1c1,n-1,a,b+1,c+1);
% derivatives and interp ops in deriv basis
Dx = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dy = D1_tri(a,b,c,H_abc,H_ab1c1,1);
% make coeffs of fref under (a,b,c) and derivatives
cfref_abc = V_abc'*FrefW;
cdxf_a1bc1 = Dx*cfref_abc;
cdyf_ab1c1 = Dy*cfref_abc;
% check representation
disp(norm(V_abc*cfref_abc-Fref)/normFref);
disp(norm(V_a1bc1*cdxf_a1bc1-df(R,S))/norm(df(R,S)));
disp(norm(V_ab1c1*cdyf_ab1c1-df(R,S))/norm(df(R,S)));

% reference tri
Rv = [0,1,0];
Sv = [0,0,1];

%img = imread("mathias-reding-vU9-VO-4Nk0-unsplash.jpg");
img = imread("cat_greyscale.jpg");
img = rgb2gray(img);
imgsz = size(img);
nXpix = imgsz(1); nYpix = imgsz(2);
[Xpix,Ypix] = meshgrid(1:nXpix,1:nYpix); 

nSamp = 20;
[X,Y] = meshgrid(linspace(1,nXpix,nSamp), linspace(1,nYpix,nSamp));
X = X(:); Y = Y(:);
DT = delaunayTriangulation(X,Y);
triInd = DT.pointLocation(Xpix(:),Ypix(:));
T = DT.ConnectivityList;
figure(1);
triplot(T,X,Y); hold on;
imagesc(img','AlphaData',0.5)
nTri = length(T);
disp(nTri);
interpolator = @(x,y) interp2(Xpix,Ypix,double(img)',x,y,'cubic');
Xpix = Xpix(:); Ypix = Ypix(:);
Imgapprox = zeros(size(img'));
ImgGradNorm = zeros(size(img'));
% TODO: Fix division by zero in jpoly_tri - cancels analytically 
% when y=1-x
for j = 1:nTri

    Xe = [X(T(j,:))';Y(T(j,:))'];
    Ixe = IncidenceMatrix(Xe);
    XYe = (Ixe * [R,S]' + Xe(:,1))';
    imginterp = interpolator(XYe(:,1),XYe(:,2));
    cimg = V_abc' * (imginterp .* W);
    cdximg = Dx*cimg;
    cdyimg = Dy*cimg;
    % find pixels inside current triangle
    xpix = Xpix(triInd==j);
    ypix = Ypix(triInd==j);
    % map those pixels to reference triangle
    XYpix = Ixe \ [xpix'-Xe(1,1);ypix'-Xe(2,1)];
    % evaluate vandermonde at pixels in reference
    Vabcpix = jPoly_tri(XYpix(1,:)',XYpix(2,:)',H_abc,n-1,a,b,c);
    Va1bc1pix = jPoly_tri(XYpix(1,:)',XYpix(2,:)',H_a1bc1,n-1,a+1,b,c+1);
    Vab1c1pix = jPoly_tri(XYpix(1,:)',XYpix(2,:)',H_ab1c1,n-1,a,b+1,c+1);
    % evaluate approximation to image within triangle
    imgapprox = Vabcpix*cimg;
    imgGradNorm = ((Va1bc1pix*cdximg).^2 + (Vab1c1pix*cdyimg).^2).^(1./2.);

    Imgapprox(triInd==j) = uint8(imgapprox);
    ImgGradNorm(triInd==j) = imgGradNorm;
    % get actual pixel values within triangle
    imgactual = interpolator(xpix,ypix);
    disp(norm(round(imgapprox)-imgactual,'fro')/norm(imgactual,'fro'));
    %plot(XYe(:,1),XYe(:,2),'k.');
    %drawnow; pause(0.01);
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
