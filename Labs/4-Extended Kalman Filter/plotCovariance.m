function plotCovariance (C , x0 , y0 , k , TYPE)

mu= [x0;y0];
NumberOfPoints = 100;
[V,D]  = eig(C);
[~,ix] = sort(diag(D), 'descend');
D      = D(ix,ix);
V      = V(:,ix);
confidence   = 2*normcdf(k) - 1;
scaling  = chi2inv(confidence, 2);
t = linspace(0, 2*pi, NumberOfPoints)';
e = [cos(t) sin(t)];
VV    = V*sqrt(D*scaling);
ee = e*(VV') + repmat(mu,1,NumberOfPoints)';

switch TYPE
    case 'Initialization'
        plot(ee(:,1), ee(:,2),'color','k');
    case 'Prior'
        plot(ee(:,1), ee(:,2),'color','r');
    case 'Update'
        plot(ee(:,1), ee(:,2),'color','g');
end

end