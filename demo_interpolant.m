function demo_interpolant
% -------------------------------------------------------------------------
% In this demo we compute the error with two test functions of the 
% interpolant construct by integrals on disks with AFS or DLS supports 
% 
% -------------------------------------------------------------------------
% Dates
% --------------------------------------------------------------------------
% First version: August 15, 2024;
% Checked: October 10, 2024.
% --------------------------------------------------------------------------
% Authors
% --------------------------------------------------------------------------
% L. Bruni Bruno and G. Elefante
% --------------------------------------------------------------------------
% Paper
% --------------------------------------------------------------------------
% "Uniform approximation of diffused data"
% L. Bruni Bruno and G. Elefante
% --------------------------------------------------------------------------
% Parameters 
degSet = [1:20];
type_extraction = 1; % 1- AFS, 2-DLS
do_plot = 1;
type_fun = 2;
ptsset = 1;
%--------------------------------------------------------------------------

M = 100;

switch ptsset 
    case 1
        x = linspace(-1,1,M);
        [xx,yy] = meshgrid(x); xx = xx(:); yy = yy(:);
        id = xx.^2 + yy.^2 <1;
        c1 = [xx(id),yy(id)];
    case 2
        M = M^2;
        c1 = rand(M,2);
        c1 = c1*2-1; 
        c1 = [c1(:,1).*sqrt(1-c1(:,2).^2/2),c1(:,2).*sqrt(1-c1(:,1).^2/2)];
end

rho = vecnorm(c1,2,2);
r1 = 1e-6.*ones(size(rho));

h = 1;
fprintf('\n\n')

for deg = degSet
    
fprintf('...................................................')
fprintf('\n \t degree          : %5d \n',deg)

dimP = nchoosek(deg+2,2);

V1 = chebVandInt(deg,c1,r1);
% V = VandInt(deg,centers,radii);

switch type_extraction
    case 1
        % AFS
        w = V1'\ones(size(V1,2),1);

        ind = find(w~=0);

        centers = c1(ind,:);
        radii = r1(ind);
    case 2
        % DLS
        [L,U,P] = lu(V1);

        idx = [1:size(V1,1)]';
        idx = P*idx;

        ind = idx(1:size(V1,2));

        centers = c1(ind,:);
        radii = r1(ind,:);
end

C = exp(1i*linspace(0,2*pi,150));

if do_plot
    figure(1)
    plot(centers(:,1),centers(:,2),'.','Color',[147,147,147]/255);
    axis equal
    hold on;
    plot(real(C),imag(C),'Color',[57,57,57]/255);
end

if do_plot
    for j=1:size(centers,1)
        CC = C*radii(j)+centers(j,1)+1i*centers(j,2);
        plot(real(CC),imag(CC),'Color',[245,119,34]/255)
    end
end

Vcheb = chebVandInt(deg,centers,radii);
% Vcheb2 = chebVandInt(deg,centers2,radii2);

switch type_fun 
    case 1
        f = @(x,y) exp(x).*sin(x+y);
    case 2
        f = @(x,y) 1./(25*(x.^2+y.^2)+1);
end

[XYW,~] = cub_unit_disk_sets(77);

for i = 1:length(radii)
    F1(i,1) = f(XYW(:,1)*radii(i)+centers(i,1),XYW(:,2)*radii(i)+centers(i,2))'*XYW(:,3)*radii(i)^2;
end

% for i = 1:length(radii2)
%     F2(i,1) = f(XYW(:,1)*radii2(i)+centers2(i,1),XYW(:,2)*radii2(i)+centers2(i,2))'*XYW(:,3)*radii2(i)^2;
% end
% F = F';

a1 = Vcheb\F1;
% a2 = Vcheb2\F2;

x = linspace(-1,1,100);
% x(1) = []; x(end) = [];
[X_old,Y_old] = meshgrid(x);
X = X_old(:).*sqrt(1-Y_old(:).^2/2); Y = Y_old(:).*sqrt(1-X_old(:).^2/2);

PolDeg = polydeg(deg); 

VxCheb = chebpolys(deg,X);
VyCheb = chebpolys(deg,Y);
for j = 1:dimP
        V(:,j) = VxCheb(:,PolDeg(j,1)+1).*VyCheb(:,PolDeg(j,2)+1);
end

P1 = V*a1;
% P2 = V*a2;

P1 = reshape(P1,size(X_old));
% P2 = reshape(P2,size(X_old));
X = reshape(X,size(X_old)); Y = reshape(Y,size(X_old));

if do_plot
    figure(4)
    mesh(X,Y,P1,'FaceAlpha','0.8','edgecolor','none','facecolor',[174,208,232]/255)
    hold on;
    mesh(X,Y,f(X,Y),'FaceAlpha','0.8','edgecolor','none','facecolor',[232,198,174]/255)
end


err1(h) = norm(P1(:)-f(X(:),Y(:)),'inf');
% err2(h) = norm(P2(:)-f(X(:),Y(:)),'inf');

fprintf('\n \t error           : %1.4e \n',err1(h))
% fprintf('\n \t error           : %1.4e \n',err2(h))
fprintf('...................................................\n')
h = h + 1;
end

if do_plot
    figure(5)
    semilogy(degSet,err1,'-o','Color',[0,142,142]/255,'LineWidth',1.5) % blue
    hold on
%     semilogy(degSet,err1,'-o','Color',[142,0,0]/255,'LineWidth',1.5) % red
%     semilogy(degSet,err,'-o','Color',[142,142,0]/255,'LineWidth',1.5) % mustard
%     hold on;
end