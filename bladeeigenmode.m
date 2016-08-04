function [ omegas, uyev, uzev ] = bladeeigenmode( radius,  twist, eixx, eiyy, mass, pitch )
% Find Rotor Blade Eigenmodes/Eigenvalues %
%  Slender beam theory by M.O. Hansen     %
% eixx = E*Ixx for element
% eiyy = E*Iyy for element
% mass = mass of element 

ni = max(size(radius));


py = zeros(ni,1);
pz = py;
flexm = zeros(2*ni-2,2*ni-2);
massmatrix = flexm;
% Flexibility matrix, y loading
    for j = 2:ni 
        m = 2*j-3;
        py(j) = 1.0;
        pz(j) = 0.0;
    [uy,~,~,uz,~,~]=deflec2(radius,eixx,eiyy,twist,pitch,py,pz);
       py(j)=0.0;
       pz(j)=0.0; 
       for i = 2:ni %Save into flex. matrix
           n = 2*i-3;
          flexm(n,m) = uy(i);
          flexm(n+1,m) = uz(i);
       end
    end
% Flexibility matrix, z loading
    for j = 2:ni 
        m = 2*j-2;
        py(j) = 0.0;
        pz(j) = 1.0;
    [uy,~,~,uz,~,~]=deflec2(radius,eixx,eiyy,twist,pitch,py,pz);
       py(j)=0.0;
       pz(j)=0.0; 
       for i = 2:ni %Save into flex. matrix
           n = 2*i-3;
          flexm(n,m) = uy(i);
          flexm(n+1,m) = uz(i);
       end
    end

% Mass Matrix %
    for i = 2:ni
        n = 2*i-3;
        massmatrix(n,n) = mass(i);
        massmatrix(n+1,n+1) = mass(i);
    end

  [V,D] = eigs(flexm*massmatrix,ni);
  %Evalues/Evectors%
    omegas = sqrt(1./(D*ones(ni,1)));
    uyev = zeros(ni,ni); uzev = uyev;
 for i=2:ni;
     n=2*i-3;
     uyev(i,1:ni)=-V(n,1:ni);
     uzev(i,1:ni)=-V(n+1,1:ni);
 end
 
 
end

