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
S = importdata("../bin/ytri_N325_n24_M946_m42.txt"); N = length(S);
R = importdata("../bin/xtri_N325_n24_M946_m42.txt");
W = importdata("../bin/wtri_N325_n24_M946_m42.txt");
%R = Zk(1:N); S = Zk(N+1:2*N); W = Zk(2*N+1:3*N);
% test function
f = @(x,y)(x+y).^21;
df = @(x,y) 21*(x+y).^20;
% eval test function on quadrature nodes
Fref = f(R,S); FrefW = Fref.*W; normFref = norm(Fref);
m = 21; n = m+1;
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

img = imread("mathias-reding-vU9-VO-4Nk0-unsplash.jpg");
%img = imread("cat_greyscale.jpg");
%img = imread('ngc6543a.jpg');
%img = imread('Parakeeets.jpg');
img = rgb2gray(img); 
nsamp = 20; nruns = 5;
[DT, triInd, XXpix, YYpix] = triangulate(img, nsamp, nruns, false);
%[DT, triInd, XXpix, YYpix] = triangulate_fast(img, nsamp, nruns);

T = DT.ConnectivityList;
X = DT.Points(:,1); Y = DT.Points(:,2);
figure(1);
triplot(T,X,Y); hold on;
imagesc(img','AlphaData',0.5)

nTri = length(T);
disp(nTri);
interpolator = @(x,y) interp2(XXpix,YYpix,double(img)',x,y,'cubic');
Xpix = XXpix(:); Ypix = YYpix(:);
Imgapprox = zeros(size(img'));
%ImgGradNorm = zeros(size(img'));
for j = 1:nTri

    Xe = [X(T(j,:))';Y(T(j,:))'];
    Ixe = IncidenceMatrix(Xe);
    XYe = (Ixe * [R,S]' + Xe(:,1))';
    imginterp = interpolator(XYe(:,1),XYe(:,2));
    cimg = V_abc' * (imginterp .* W);
    %cdximg = Dx*cimg;
    %cdyimg = Dy*cimg;
    % find pixels inside current triangle
    xpix = Xpix(triInd==j);
    ypix = Ypix(triInd==j);
    % map those pixels to reference triangle
    XYpix = Ixe \ [xpix'-Xe(1,1);ypix'-Xe(2,1)];
    % evaluate vandermonde at pixels in reference
    Vabcpix = jPoly_tri(XYpix(1,:)',XYpix(2,:)',H_abc,n-1,a,b,c);
    %Va1bc1pix = jPoly_tri(XYpix(1,:)',XYpix(2,:)',H_a1bc1,n-1,a+1,b,c+1);
    %Vab1c1pix = jPoly_tri(XYpix(1,:)',XYpix(2,:)',H_ab1c1,n-1,a,b+1,c+1);
    
    % evaluate approximation to image within reference triangle
    imgapprox = Vabcpix*cimg;
    %imgGradNorm = ((Va1bc1pix*cdximg).^2 + (Vab1c1pix*cdyimg).^2).^(1./2.);

    Imgapprox(triInd==j) = uint8(imgapprox);
    %ImgGradNorm(triInd==j) = imgGradNorm;
    % get actual pixel values within triangle
    %imgactual = interpolator(xpix,ypix);
    %err = norm(imgapprox-imgactual,'fro')/norm(imgactual,'fro');
    %plot(XYe(:,1),XYe(:,2),'k.');
    %drawnow; pause(0.01);
    %disp([j,err]);img
    disp(j);
end


function Ixe = IncidenceMatrix(Xe)

Ixe = [Xe(:,2)-Xe(:,1), Xe(:,3)-Xe(:,1)];

end

function [DT, triInd, XXpix, YYpix] = triangulate(I, nSamp, nRuns, go)

% first insert evenly spaced candidate coordinates
% from subsampled version of I
imgsz = size(I);
nXpix = imgsz(1); nYpix = imgsz(2);
[XXpix,YYpix] = meshgrid(1:nXpix,1:nYpix); 
interpolator = @(x,y) interp2(XXpix,YYpix,double(I)',x,y,'linear');
Xpix = XXpix(:); Ypix = YYpix(:);
[X,Y] = meshgrid(linspace(1,nXpix,nSamp), linspace(1,nYpix,nSamp));
% now create Delaunay triangulation
X = X(:); Y = Y(:);
% load low order quad rule
R = importdata("../bin/xtri_N55_n9_M91_m12.txt");
S = importdata("../bin/ytri_N55_n9_M91_m12.txt");
W = importdata("../bin/wtri_N55_n9_M91_m12.txt");
m = 6; n = m+1;
a = 0.5; b = a; c = a;
H_abc = structure_factors_tri(n+1,a,b,c);
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
for iRun = 1:nRuns
  DT = delaunayTriangulation(X,Y);
  % triInd is list of triangle indices
  % for each point (Xpix,Ypix), corresponding
  % to that point exisiting in triangle
  triInd = DT.pointLocation(Xpix,Ypix);
  % now we generate I'
  T = DT.ConnectivityList;
  triplot(T,X,Y); drawnow; pause(0.5);
  if ~go
    break;
  end
  nTri = length(T);

  nP = length(DT.Points);
  PSNRS = zeros(nP,1);
  % proceed vertex by vertex
  for ip = 1:nP
    % remove the vertex if not on boundary
    Xp = X; 
    Yp = Y; 
    if (Xp(ip) ~= 1 && Xp(ip) ~= nXpix && ...
        Yp(ip) ~= 1 && Yp(ip) ~= nYpix)
      Xp(ip) = [];
      Yp(ip) = [];
      DTp = delaunayTriangulation(Xp,Yp);
      Tp = DTp.ConnectivityList;
      nTrip = length(Tp);
      triIndp = DTp.pointLocation(Xpix(:),Ypix(:));
      Ipp = zeros(size(I'));
      for j = 1:nTrip
        % Get vertices of current triangle
        Xe = [Xp(Tp(j,:))';Yp(Tp(j,:))'];
        Ixe = IncidenceMatrix(Xe);
        % map reference quad points to current triangle
        XYe = (Ixe * [R,S]' + Xe(:,1))';
        % linearly interpolate to quad points in current triangle
        Iinterp = interpolator(XYe(:,1),XYe(:,2));
        % compute coefficients
        cimg = V_abc' * (Iinterp .* W);
        % find pixels inside current triangle
        xpix = Xpix(triIndp==j);
        ypix = Ypix(triIndp==j);
        % map those pixels to reference triangle
        XYpix = Ixe \ [xpix'-Xe(1,1);ypix'-Xe(2,1)];
        % evaluate vandermonde at pixels in reference
        Vabcpix = jPoly_tri(XYpix(1,:)',XYpix(2,:)',H_abc,n-1,a,b,c);
        % evaluate approximation to image within reference triangle
        imp = Vabcpix*cimg;
        Ipp(triIndp==j) = imp;
      end
      % compute PSNR between original I and approximate Ip
      PSNRp = psnr(uint8(Ipp),I');
      PSNRS(ip) = PSNRp;
      disp([ip, PSNRp]);
    end
  end
  [~,maxInd] = max(PSNRS);
  X(maxInd) = []; Y(maxInd) = [];
  disp(iRun);
end

end

function [DT, triInd, XXpix, YYpix] = triangulate_fast(I, nSamp, nRuns)

% first insert evenly spaced candidate coordinates
% from subsampled version of I
imgsz = size(I);
nXpix = imgsz(1); nYpix = imgsz(2);
[XXpix,YYpix] = meshgrid(1:nXpix,1:nYpix); 
interpolate = @(x,y) interp2(XXpix,YYpix,double(I)',x,y,'nearest');
Xpix = XXpix(:); Ypix = YYpix(:);
[X,Y] = meshgrid(linspace(1,nXpix,nSamp), linspace(1,nYpix,nSamp));
X = X(:); Y = Y(:); 
for iRun = 1:nRuns
  nP = length(X);
  PSNRS = zeros(nP,1);
  % proceed vertex by vertex
  for ip = 1:nP
    % remove the vertex if not on boundary
    Xp = X;
    Yp = Y; 
    if (Xp(ip) ~= 1 && Xp(ip) ~= nXpix && ...
        Yp(ip) ~= 1 && Yp(ip) ~= nYpix)
      Xp(ip) = [];
      Yp(ip) = [];
      % compute temp Delaunay triangulation
      DTp = delaunayTriangulation(Xp,Yp);
      Tp = DTp.ConnectivityList;
      nTrip = length(Tp);
      triIndp = DTp.pointLocation(Xpix(:),Ypix(:));
      Ipp = zeros(size(I'));
      % approximate image
      for j = 1:nTrip
        % Get vertices of current triangle
        Xe = [Xp(Tp(j,:))';Yp(Tp(j,:))'];
        % get color at each vertex of current triangle
        vcols = interpolate(Xe(1,:)',Xe(2,:)');
        % find pixels inside current triangle
        xpix = Xpix(triIndp==j);
        ypix = Ypix(triIndp==j);
        % get barycentric weights for those pixels
        [Wv1,Wv2,Wv3] = baryweights(Xe,xpix,ypix);
        % get  colors within triangle via barycentric interp
        imp = vcols(1)*Wv1 + vcols(2)*Wv2 + vcols(3)*Wv3;
        Ipp(triIndp==j) = imp;
      end
      % compute PSNR between original I and approximate Ip
      PSNRp = ssim(uint8(Ipp),I');
      PSNRS(ip) = PSNRp;
      disp([ip, PSNRp]);
    end
  end
  [~,maxInd] = max(PSNRS);
  X(maxInd) = []; Y(maxInd) = [];
  DT = delaunayTriangulation(X,Y);
  % triInd is list of triangle indices
  % for each point (Xpix,Ypix), corresponding
  % to that point exisiting in triangle
  triInd = DT.pointLocation(Xpix,Ypix);
  % now we generate I'
  triplot(DT.ConnectivityList,X,Y); drawnow; 
  disp(iRun);
end

end

function [Wv1,Wv2,Wv3] = baryweights(Xe,X,Y)

Xv1 = Xe(1,1); Xv2 = Xe(1,2); Xv3 = Xe(1,3);
Yv1 = Xe(2,1); Yv2 = Xe(2,2); Yv3 = Xe(2,3);
denom = (Yv2-Yv3)*(Xv1-Xv3) + (Xv3-Xv2)*(Yv1-Yv3);
Wv1 = ((Yv2-Yv3)*(X-Xv3) + (Xv3-Xv2)*(Y-Yv3))/denom;
Wv2 = ((Yv3-Yv1)*(X-Xv3) + (Xv1-Xv3)*(Y-Yv3))/denom;
Wv3 = 1-Wv1-Wv2;

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


%   Ip = zeros(size(I'));
% 
%   for j = 1:nTri
%     % Get vertices of current triangle
%     Xe = [X(T(j,:))';Y(T(j,:))'];
%     Ixe = IncidenceMatrix(Xe);
%     % map reference quad points to current triangle
%     XYe = (Ixe * [R,S]' + Xe(:,1))';
%     % linearly interpolate to quad points in current triangle
%     Iinterp = interpolator(XYe(:,1),XYe(:,2));
%     % compute coefficients
%     cimg = V_abc' * (Iinterp .* W);
%     % find pixels inside current triangle
%     xpix = Xpix(triInd==j);
%     ypix = Ypix(triInd==j);
%     % map those pixels to reference triangle
%     XYpix = Ixe \ [xpix'-Xe(1,1);ypix'-Xe(2,1)];
%     % evaluate vandermonde at pixels in reference
%     Vabcpix = jPoly_tri(XYpix(1,:)',XYpix(2,:)',H_abc,n-1,a,b,c);
%     % evaluate approximation to image within reference triangle
%     ip = Vabcpix*cimg;
%     Ip(triInd==j) = ip;
%   end
% 
%   % compute PSNR between original I and approximate Ip
%   PSNR = psnr(I',uint8(Ip));
%   disp(PSNR);
% 
