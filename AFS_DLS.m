clear all

ptsset = 1;
degV = [1:15];

M = 100;

switch ptsset 
    case 1
        x = linspace(-1,1,M);
        [xx,yy] = meshgrid(x); xx = xx(:); yy = yy(:);
        id = xx.^2 + yy.^2 <1;
        centers = [xx(id),yy(id)];
    case 2
        M = M^2;
        centers = rand(M,2);
        centers = centers*2-1; 
        centers = [centers(:,1).*sqrt(1-centers(:,2).^2/2),centers(:,2).*sqrt(1-centers(:,1).^2/2)];
end

rho = vecnorm(centers,2,2);
radii = 1e-6.*ones(size(rho));

hh = 1;
l1 = zeros(size(degV,2),1); l2 = zeros(size(degV,2),1);
for deg = degV

    fprintf('\n...................................................')
    fprintf('\n \t degree          : %5d \n',deg)
    V = chebVandInt(deg,centers,radii);
    % V = VandInt(deg,centers,radii);
    dim(hh) = nchoosek(deg+2,2);


    % AFS
    w = V'\ones(size(V,2),1);

    ind = find(w~=0);

    c1 = centers(ind,:);
    r1 = radii(ind);

    l1(hh) = LebesgueConstant(deg,c1,r1);

    % DLS
    [L,U,P] = lu(V);

    idx = [1:size(V,1)]';
    idx = P*idx;

    ind = idx(1:size(V,2));

    c2 = centers(ind,:);
    r2 = radii(ind,:);

    l2(hh) = LebesgueConstant(deg,c2,r2);


    fprintf('\n...................................................')
    fprintf('\n \t # AFS  : %5d   \t # DLS  : %5d',size(c1,1), size(c2,1))
    fprintf('\n \t Leb.C. : %3.2e \t Leb.C. : %3.2e ',l1(hh),l2(hh))
    fprintf('\n...................................................')
    hh = hh+1;
end
fprintf('\n')

figure(1)
semilogy(degV,l1,'-o','Color',[142,0,0]/255,'LineWidth',1.5)
hold on
semilogy(degV,dim,'--','Color',[142,142,142]/255,'LineWidth',1.5)
xlim([degV(1),degV(end)])
ylim([1 dim(end)])

figure(2)
semilogy(degV,l2,'-o','Color',[142,0,0]/255,'LineWidth',1.5)
hold on
semilogy(degV,dim,'--','Color',[142,142,142]/255,'LineWidth',1.5)
xlim([degV(1),degV(end)])
ylim([1 dim(end)])