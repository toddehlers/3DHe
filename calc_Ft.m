function [Fts,V,S,Rs] = calc_Ft(mineral,geometry,param)

% function to calculate Ft-values according to Ketcham et al. (2011): 
% Accounting for long alpha-particle stopping distances
% in (U?Th?Sm)/He geochronology: Refinement of the baseline case

% load mean stopping distances (SRIM-2003)
load stopping_distances.mat



% possible geometries: 

% ellipsoid
% required measurements:    a,b,c = semi-principal axes (a<b<c)

% cylinder
% required measurements:    h = height 
%                           r = radius

% block
% required measurements:    a,b,c = principal axes (a<b<c)
%                           Np = number of pyramids (0,1,2)

% hexagonal
% required measurements:    H = complete length along c-axis
%                           W = maximum width (from edge to edge)
%                           L = minimum width (from face to face)
%                           Np = number of pyramids (0,1,2)

% triclinic
% required measurements:    a,b,c = principal axes (a<b<c)
%                           alpha,beta,gamma = angles between (b,c), (c,a), (a,b)

% monoclinic
% required measurements:    c = complete height along c-axis (with pyramids)
%                           b,c = principal axes (a<b)
%                           gamma = angle between (b,c)
%                           Np = number of pyramids (1 or 2)

% transfer variables
switch geometry
    case 'ellipsoid (A)'
        a=str2num(param.a{1});
        b=str2num(param.b{1});
        c=str2num(param.c{1});
    case 'cylinder (B)'
        h=str2num(param.h{1});
        r=str2num(param.r{1});
    case 'block (C)'
        c=str2num(param.a{1});
        b=str2num(param.b{1});
        a=str2num(param.c{1});
        Np=param.Np;
    case 'hexagonal (D)'
        H=str2num(param.H{1});
        W=str2num(param.W{1});
        L=str2num(param.L{1});
        Np=param.Np;
    case 'triclinic (E)'
        a=str2num(param.a{1});
        b=str2num(param.b{1});
        c=str2num(param.c{1});
        alpha=str2num(param.alpha{1});
        beta=str2num(param.beta{1});
        gamma=str2num(param.gamma{1});
    case 'monoclinic (F)'
        a=str2num(param.a{1});
        b=str2num(param.b{1});
        c=str2num(param.c{1});
        gamma=str2num(param.gamma{1});
        Np=param.Np;
end

switch mineral
    case 'ap'
        stop_dist=apatite;
    case 'zr'
        stop_dist=zircon;
end
isotops=fieldnames(stop_dist);
for i=1:length(isotops)
    % R is the mean stoping distance of 238, 235U, 232Th or 147 Sm in apatite or zircon
    R=stop_dist.(isotops{i});
    switch geometry
        case {'ellipsoid','ellipsoid (A)'}
            if a/b<0.3 || a/b>1 || c/b<1 || c/b>10
                warning('Analytical function not applicable, Ft may be inaccurate!') 
            end
            if a==b==c
                Rs=a;
            else
                V=4/3*pi()*a*b*c;
                % surface approximation with Knud Thomsen's formula
                p=1.6075;
                S=4*pi()*((a^p*b^p+b^p*c^p+c^p*a^p)/3)^(1/p);
                Rs=3*V/S;
            end
            Ft=1-3/4*R/Rs+(1/16+0.1686*(1-a/Rs)^2)*(R/Rs)^3;
            %Ft=1-3/4*R/Rs+(1/16+0.1686*(1-a/Rs)^2)*(1.31*R^3)/Rs^3;
        case {'cylinder','cylinder (B)'}
            if h/r<0.2 || h/r>2000
                 warning('Analytical function not applicable, Ft may be inaccurate!') 
            end
            Rs=3/2*r*h/(r+h);
            Ft=1-1/2*((r+h)*R)/(r*h)+0.2122*(1.09*R^2/(r*h))+0.0153*1.31*R^3/r^3;
        case {'block','block (C)'}
            if a<2*R || b<2*R || c<2*R || a/b<0.3 || a/b>1 || c/b<1 || c/b>10
                 warning('Analytical function not applicable, Ft may be inaccurate!') 
            end
            if Np==0
                Rs=3/2*(1/(1/a+1/b+1/c));
                Ft=1-3/4*R/Rs+0.2122*(a+b+c)*1.09*R^2/(a*b*c)-0.00995*1.09*R^3/(a*b*c);
            else
                V=a*b*c-Np*min([a b])/4*(max([a b])^2+min([a b])^2/3);
                S=2*(a*b+b*c+a*c)-Np*((a^2+b^2)/2+(2-sqrt(2))*a*b);
                Rs=3*V/S;
                Ft=1-3/4*R/Rs+(0.2095*(a+b+c)-(0.096-0.013*(a^2+b^2)/c^2)*(a+b)*Np)*1.09*R^2/V;
            end
        case {'hexagonal','hexagonal (D)'}
            if H/W<0.1 || H/W>2000 || W/L>3.85 || W/L<1/sqrt(3)
                 warning('Analytical function not applicable, Ft may be inaccurate!') 
            end
            if L>sqrt(3)/2*W
                dV=1/(6*sqrt(3))*(L-sqrt(3)/2*W)^3;
            else
                dV=0;
            end
            V=H*L*(W-L/(2*sqrt(3)))-Np*(sqrt(3)/8*W^2*L-dV);
            S=2*H*(W+L/sqrt(3))+2*L*(W-L/(2*sqrt(3)))-Np*(sqrt(3)/4*W^2+(2-sqrt(2))*W*L+(sqrt(2)-1)/(2*sqrt(3))*L^2);
            Rs=3*V/S;
            Ft=1-3/4*R/Rs+((0.2093-0.0465*Np)*(W+L/sqrt(3))+(0.1062+(0.2234*R)/(R+6*(W*sqrt(3)-L)))*(H-Np*(W*sqrt(3)/2+L)/4))*1.09*R^2/V;
        case {'triclinic','triclinic (E)'}
            if a/b<0.3 || a/b>1 || c/b<1 || c/b>10
                 warning('Analytical function not applicable, Ft may be inaccurate!') 
            end
            d=1-cosd(alpha)^2-cosd(beta)^2-cosd(gamma)^2+2*cosd(alpha)*cosd(beta)*cosd(gamma);
            V=a*b*c*sqrt(d);
            if d>0
                S=2*(b*c*sind(alpha)+a*c*sind(beta)+a*b*sind(gamma));
            else
                warning('Analytical function not applicable, Ft may be inaccurate!') 
            end
            Rs=3*V/S;
            Ft=1-3/4*R/Rs+(-0.1651*(a+b+c)+0.3746*(a*sind(beta)*sind(gamma)+b*sind(alpha)*sind(gamma)+c*sind(alpha)*sind(beta))/sqrt(d))*1.09*R^2/(a*b*c*sqrt(d));
        case {'monoclinic','monoclinic (F)'}
            V=a*b*c*sind(gamma)-Np*min([a b])/4*(max([a b])^2+min([a b])^2/3)*sind(gamma)^2;
            S=2*(a*c+b*c+a*b*sind(gamma))-Np*((a^2+b^2)/2+(2-sqrt(2))*a*b)*sind(gamma);
            Rs=3*V/S;
            Ft=1-3/4*R/Rs+(0.2095*(a+b+c)+0.3746*(1/sind(gamma)-1)*c-(0.051+0.045*sind(gamma)+(0.069-0.082*sind(gamma))*(a^2+b^2)/c^2)*(a+b)*Np)*1.09*R^2/V;
        otherwise
            warning('Geometry not known!') 
    end
    Fts.(isotops{i})=Ft;
end
        