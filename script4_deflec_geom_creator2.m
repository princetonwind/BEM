% Script for Finding Blade Deflection %
% Solid blade material only           %
clear all
close all
% Rotor Blade Inputs %
    E = 200 * 10^9;     %Young's modulus (Pa)
    rho = 7950;         %Blade density, kg/m^3

% Rotor geometry file location %
    rtype = 'C:\Users\tnealon\Documents\Princeton Research\V27 Data\V27_Full-Scale_BEM_Geometry.txt';
    fid = fopen(rtype,'r');
    rotorgem = textscan(fid,'%f %f %f %s','Delimiter',' ','MultipleDelimsAsOne',1,'Headerlines',1);
    fclose(fid);
    rotor = [rotorgem{1} rotorgem{2} rotorgem{3}];
    foil = rotorgem{4};
    
%% Approximate Area Moments of Intertia and Mass %%
     Ixx = zeros(numel(rotor(:,1)),1); Iyy = Ixx; area = Ixx; 
     ymin = Ixx; ymax = Ixx;
for i = 1:numel(rotor(:,1))
    disp(['Current Foil: ' char(foil(i))])
   [Ixx(i),Iyy(i), Atot] = ainertia(char(foil(i)),rotor(i,2)); 
    %Find mass and distribution%
    area(i) = Atot;
end

mass = zeros(numel(rotor(:,1))-1,1); rcen = mass; rav = mass;
%Find mass and centroids%
for n = 1:numel(rotor(:,1))-1
    r2 = rotor(n+1,1); r1 = rotor(n,1);
   phi = (area(n+1) - area(n))/(r2 - r1);
   mass(n) = rho * (area(n)*(r2 - r1) + ...
       phi / 2 * (r2^2 + r1^2) - phi * r2 * r1);
   % center of mass for each element, csys at beginning of element     
   rcen(n) = rho / mass(n) * (area(n)/2 * (r2^2 - r1^2) ...
       - phi/2 * r1 * r2^2 + phi/3 * (r2^3 + r1^3/2));
%    rcen(n) = rcen(n) + rotor(n); %shift to blade coordinate system
   rav(n) = (r1 + r2) / 2;
end
   

figure(3)
plot(rcen,mass,'rx-',rav,mass,'bo-',rotor(:,1),ones(numel(rotor(:,1)),1),...
    'gs','Markersize',8)
legend('Centroid','average distance','Radial AFs')

    [wfolder,name,~] = fileparts(rtype);
    fopt = [wfolder '\' name '_deflec_geom.txt'];
        radius = rcen; 
        %Assume linear variation for all properties
        chord = interp1(rotor(:,1),rotor(:,2),rcen); 
        twist = interp1(rotor(:,1),rotor(:,3),rcen);
        EIxx  = interp1(rotor(:,1),E.*Ixx,rcen); 
        EIyy  = interp1(rotor(:,1),E.*Iyy,rcen);
%         Area  = interp1(rotor(:,1),area,rcen);
        T = table(radius,chord,twist,EIxx,EIyy,mass);
    writetable(T,fopt,'FileType','text','delimiter','\t')

    disp(T)

  figure(2)
  subplot(2,1,1)
  plot(radius,EIxx,'b--',radius,EIyy,'k--',...
      'Linewidth',2)
  ylabel('E*I (m^4)')
  legend('EIxx','EIyy')
  subplot(2,1,2)
  plot(radius,mass.*1000,'r-','Linewidth',2)
  xlabel('Radius (m)')
  ylabel('Mass (grams)')
  
  
  
  
    
    