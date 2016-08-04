function [ rho, mu, summ, bldout ] = mark_bem360(rotor, foil, T, P, speed, pitch , U , nb, NE)
%Blade Element Momentum code for HRTF
% rotor = [radius (m), chord (m), twist (deg)]
% ROTOR and FOIL must be in order of increasing radius!!
% foil = { 'cell of airfoil strings'  }
% T = Tunnel Temp (degC); P = Ambient (Tunnel) Pressure (Pa)
% speed = Rotation Speed (RPM);   pitch = blade pitch (deg)
% nb = # of blades; NE = # Blade Elements; 
% relax = relaxation factor as in Qblade documentation
% 
%%% Mark Miller 11-25-15 %%% :]
%   

% location of the airfoil data 
rfolder = 'C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data';
nmax = 1000; %Maximum number of iterations allowed
tol = 1E-6;
relax = 0.5;

     %  Check to make sure input values are all the correct size  %
      if numel(T) ~= numel(P) || (numel(T) ~= numel(speed) || numel(speed) ~= numel(P))
          error('Input Conditons not equal size (T,P,speed)!!')
      elseif numel(U) > 1
          error('Only one free-stream velocity value can be run!')
      elseif numel(foil) ~= numel(rotor(:,1))
          disp(['Number of airfoils: ' num2str(numel(foil))])
          disp(['Number of radial locations: ' num2str(numel(rotor(:,1)))])
          error('Airfoils not assigned or too few radius locations given')
      end
    
TSR = speed*pi.*rotor(end,1)./(30*U);   
   rho = 1.225;
   mu = 1.6*10^-5;
% Discretize blade in Blade Elements %
    % Find length of blade elements
    db = (rotor(end,1) - rotor(1,1)) / NE;
    % create array of radial locations for blade elements  
    abr = linspace(rotor(1,1) + db/2, rotor(end,1) - db/2, NE)'; 
    % Interpolate array of chord elements
    abc = interp1(rotor(:,1),rotor(:,2),abr);
    % Interpolate twist for elements, in degrees
    abt = interp1(rotor(:,1),rotor(:,3),abr);
    % Find airfoils for each blade element
    Afidx = round(interp1(rotor(:,1),1:1:length(rotor(:,1)-1),abr));
    for mm = 1:numel(Afidx)
        %populate list with airfoils
       aft(mm) = foil(Afidx(mm)); 
    end    

%  Determine tunnel Conditions  %
      %[rho, mu ] = Comp_air_corr(T,P);
      % Find local pitch angle in radians %
      theta = (abt + pitch).*pi./180;
      % Find rotor Solidity %
      sigma = abc * nb ./ (2*pi .* abr);

% Load airfoil data  
        for m = 1:numel(aft)
            %Load data into 3-dimensional matrix, 
            % ( AoA , Re , Foil Type )
            dum = dlmread([rfolder '\' char(aft(m)) '_CL.txt'],'\t',1,0);
            Cl(:,:,m) = dum(2:end,2:end);
            Re(:,m) = dum(1,2:end);
            aoa(:,m) = dum(2:end,1).*pi./180; %Convert to radians
            Cd(:,:,m) = dlmread([rfolder '\' char(aft(m)) '_CD.txt'],'\t',2,1);
            clear dum
        end
%-------------------%    
%  Enter BEM Code   %
%-------------------%

      %  Reserve variables  %
      a = zeros(numel(abr),1); ap = a; alpha = a; phi = a;
      Clc = a; Cdc = a; Cn = a; Ct = a; Ulocal = a; Rec = a;

    for jj = 1:length(abr) %index through all the foils separately
        m = 1; %re-initialize while loop counter
	
        while m <= nmax + 1
            if m == nmax %Error message if maximum iteration count is reached
                disp('Maximum # iterations reached for airfoil at ')
                disp(['R = ' num2str(abr(jj)) ' and foil: ' char(aft(jj))])
                break
            end

         %Find local Reynolds # for Cl/Cd lookup
             Ulocal(jj) = ((U.*(1-a(jj))).^2 + (speed.*2*pi.*abr(jj)./60.*(1+ap(jj))).^2).^(0.5);
             Rec(jj) = rho.*abc(jj).*Ulocal(jj)./mu;
%              Rec(jj) = 0.2E6;  %Forces a specific Re
         % Find inflow angle for each position  
             phi(jj) = atan((1-a(jj)).*U./((1+ap(jj)).*speed.*pi.*abr(jj)./30));
         % Find local aoa
             alpha(jj) = phi(jj) - theta(jj);

         % Find closest Re for each section %
            dmin = abs(Re(:,jj)-repmat(Rec(jj)',[size(Re,1) 1])); 
            [Y, idx] = min(dmin);
            clear dmin   
         % Interpolate between available AoA
             Clc(jj) = interp1(aoa(:,jj),Cl(:,idx,jj),alpha(jj));
             Cdc(jj) = interp1(aoa(:,jj),Cd(:,idx,jj),alpha(jj));
         % Transform to Normal and Tangential Directions
            Cn(jj) = Clc(jj)*cos(phi(jj)) + Cdc(jj)*sin(phi(jj));
            Ct(jj) = Clc(jj)*sin(phi(jj)) - Cdc(jj)*cos(phi(jj));
         % Prandtl's Tip Loss Factor, set F=1 to turn off
            if sin(phi(jj)) < 0.02 %Avoid sigularities in correction
                F = 1;
            else
                eff = nb/2*(rotor(end,1) - rotor(jj))/(rotor(jj)*sin(phi(jj)));
                F = 2*acos(exp(-eff)) / pi;
            end 
          %Use Glauert Correction
            if a(jj) > 0.2 %From Spera (1994), see Hansen CH 6
                K = 4.*F.*sin(phi(jj))^2. / (sigma(jj) * Cn(jj));
                anew = 0.5*(2 + K*(1-2*0.2) - sqrt((K*(1-2*0.2)+2)^2 + ...
                    4*(K*0.2^2 - 1)));
            else
            anew = 1./(4*F*sin(phi(jj))^2/(sigma(jj)*Cn(jj)) + 1);
            end
            apnew = 1./(4*F*sin(phi(jj))*cos(phi(jj))/(sigma(jj)*Ct(jj)) - 1 );
          % Relaxation constant set to 1 for first few iterations as in Qblade docs    
            if m < 2 
                arelax = 1;
            elseif m == 2
                % Use 3 point eqn to find center of initial oscillations
                % See Qblade documentation
                arelax = 1;
                anew = 0.25*anew + 0.5*a(jj) +0.25*aprev; 
            else
                arelax = relax;
            end
        % Calculate error for each foil %
            if abs(anew - a(jj)) < tol 
                a(jj) =  anew;
                ap(jj) = apnew;
                break
            else
                aprev = a(jj); %Save previous iteration value
                a(jj) = arelax*anew +(1-arelax)*a(jj); 
                ap(jj) = arelax*apnew + (1-arelax)*ap(jj);
            end

        clear anew apnew idx Y
        
        m = m + 1; %Iteration Counter  
        end
    end
        %Find normal and tangential forces per unit length
            Pn = 0.5*rho.*Ulocal.^2.*abc.*Cn;
            Pt = 0.5*rho.*Ulocal.^2.*abc.*Ct;
        %Calculate hub moment, linear variation btwn points
        A = zeros(numel(abr)-1,1); B = A; M = A;
        for m = 1:numel(abr)-1
            A(m) = (Pt(m+1) - Pt(m))/(abr(m+1) - abr(m));
            B(m) = (Pt(m)*abr(m+1) - Pt(m+1)*abr(m)) / ...
                (abr(m+1) - abr(m));
            M(m) = (1/3)*A(m)*(abr(m+1)^3 - abr(m)^3) + ...
                (1/2)*B(m)*(abr(m+1)^2 - abr(m)^2);
        end
        %Calculate Thrust force loading
            Tn = trapz(abr,Pn);
            Mtot = nb * sum(M); %Total hub moment
            Ttot = nb * Tn;     %Total axial thrust
        
        %Calculate global properties:
        CT = Ttot./(0.5*rho*U^2*pi*abr(end)^2);
        CP = (Mtot*speed.*pi./30)./(0.5*rho*U^3*pi*abr(end)^2);
% whos
ReD = rho*U*rotor(end,1)*2./mu;
summ = [U TSR ReD Rec(end) CT CP Ttot Mtot*speed.*pi./30 ];
bldout = [abr Pt Pn a Clc Cdc alpha phi speed*pi/30.*abr.*(1+ap) U.*(1-a) ];
end

