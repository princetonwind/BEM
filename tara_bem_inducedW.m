%Model V27 BEM Code V 2.0
%Tara Nealon

clc, clear;

T = 20; %temperature, in C
P = 103100; %tunnel pressure
speed = 43; %rotational speed (1)
pitch = 0; %initial pitch angle
%U = 5;
U = 3.5:0.5:25; %Starting freestream velocity
Bl = 3; %number of blades
NE = 9; %number of geometry points
scale = 1.0; %scale to change to the model size- not input yet

%Bending Moment without centrifugal loading
mbend_noc = [998925.897863101,871750.108592105,691468.185178063,453372.424868382,251253.097358051,76234.8277348876,-1.09902709751445,0];

%Open Rotor Geometry file- change to V27_Model_BEM_Geometry.txt for model
file = 'C:\Users\tnealon\Documents\Princeton Research\V27 Data\V27_Full-Scale_BEM_Geometry.txt';
f_id = fopen(file,'r');
rotorgeom = textscan(f_id,'%f %f %f %s','Delimiter',' ','MultipleDelimsAsOne',1,'Headerlines',1);
fclose(f_id);
rotor = [rotorgeom{1} rotorgeom{2} rotorgeom{3}]; %Radius, Chord, Twist
foil = rotorgeom{4}; %Airfoil

folder = 'C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data';
max_iter = 1000; %Maximum number of iterations allowed
tlr = 1E-6;
relax = 0.5;

%Read EI and mass data
etype = 'C:\Users\tnealon\Documents\Princeton Research\V27 Data\V27_Full-Scale_BEM_Geometry_deflec_geom.txt';
fid = fopen(etype,'r');
eim = textscan(fid,'%f %f %f %f %f %f','Delimiter',' ','MultipleDelimsAsOne',1,'Headerlines',1);
fclose(fid);
EIxx = eim{4};
EIyy = eim{5};
mass = eim{6};

p_y = zeros(NE,1); p_z = p_y; F_x = p_y;

radius = eim{1};
twist = eim{3};
eixx = EIxx;
eiyy = EIyy;
[ omegas, uyev, uzev, thetaz, mbendy, mbendz ] = tara_bladeeigenmode(radius,  twist, eixx, eiyy, mass, pitch, p_y, p_z);

storedy = mbendy;
storedz = mbendz;
iter = 1;
%Begin BEM
p_cent = zeros(NE,1);
t = 0;
A = zeros(NE-1,1); B = A; M = A;

%Start loop for changing velocity 
for i = 1:numel(U);
    disp(['Current Tunnel Velocity ' num2str(U(i))])
    TSR = speed*pi.*rotor(end,1)./(30*U(i)); %Calculate current TSR
    disp(['Current TSR ' num2str(TSR)])
        
    %Set up matrices for read input
    %blade element length
    db = (rotor(end,1) - rotor(1,1))/NE;
    %create array of blade elements
    abe = linspace(rotor(1,1) + db/2, rotor(end,1) - db/2, NE)';
    %create array of chord for elements
    abc  = interp1(rotor(:,1),rotor(:,2),abe);
    %create array of twist for elements
    abt  = interp1(rotor(:,1),rotor(:,3),abe);
    %determine & assign airfoils for each element
    airf_id = round(interp1(rotor(:,1),1:1:length(rotor(:,1)),abe)); %got rid of the -1
    for m = 1:numel(airf_id)
        afn(m) = foil(airf_id(m));
    end
    
    %Calculate Tunnel Conditions
    [rho, mu ] = Comp_air_corr(T,P);
    
    %Calculate local pitch (blade pitch + local twist) in radians
    theta = (pitch + abt).*pi/180;
    
    %Calculate Solidity
    sigma = abc*Bl./(2*pi.*abe);
    
    %Read Airfoil Cd & Cl data
    for m = 1:numel(afn)
        %Load data into 3-dimensional matrix,
        % ( AoA , Re , Foil Type )
        dum = dlmread([folder '\' char(afn(m)) '_CL.txt'],'\t',1,0);
        Cl(:,:,m) = dum(2:end,2:end);
        Re(:,m) = dum(1,2:end);
        aoa(:,m) = dum(2:end,1).*pi./180; %Convert to radians
        Cd(:,:,m) = dlmread([folder '\' char(afn(m)) '_CD.txt'],'\t',2,1);
        clear dum
    end
    
    %Initialize variables
    a = zeros(numel(abe),1);
    tstep = 0.05;
    tmax = 30;
    a_p = a; alpha = a; phi = a; Cl_c = a; Cd_c = a; C_t = a; C_n = a; U_loc = a; Re_c = a; W_y = zeros(NE,tmax/tstep); W_z = W_y;
    W_intz = W_y; W_inty = W_y; 

    
%Time Loop- need for induced velocities 
while t < tmax
    
    %Loop through each section
    for j = 1:numel(abe)
        
        m = 1; %Reinitializes max iteration counter
    calc_er  = 10; 
       
        while calc_er >= 0.05 %m <= max_iter
            %Find local & relative velocities
            if iter == 1 %induced velocity not calculated yet
                U_loc(j,iter) = ((U(i).*(1-a(j)) ).^2 + (speed.*pi.*abe(j)./30.*(1+a_p(j))  ).^2).^(0.5); %Equivalent to V_rel_mag
                V_rel_y1(j) = (1+a_p(j)).*speed.*pi.*abe(j)./30 + W_y(j,iter);
                V_rel_z1(j) = (1-a(j)).*U(i) + W_z(j,iter);
            else
                U_loc(j,iter) = ((U(i).*(1-a(j)) + W_z(j,iter-1) ).^2 + (speed.*pi.*abe(j)./30.*(1+a_p(j)) + W_y(j,iter-1) ).^2).^(0.5); %Equivalent to V_rel_mag
                V_rel_y1(j) = (1+a_p(j)).*speed.*pi.*abe(j)./30 + W_y(j,iter-1);
                V_rel_z1(j) = (1-a(j)).*U(i) + W_z(j,iter-1);
            end
            %Calculate local Re, based on chord
            Re_c(j) = rho.*abc(j).*U_loc(j,iter)./mu;
            %Find phi (inflow angle)
            phi(j) = atan((1-a(j)).*U(i)./((1+a_p(j)).*speed.*pi.*abe(j)./30));
            %Find local AoA, in radians
            alpha(j) = phi(j) - theta (j);

            %Determine Re associated with each section, used to look up
            %Cl/Cd
            diff = abs(Re(:,j)-repmat(Re_c(j)',[size(Re,1) 1])); %Compute difference between Re and Re_c to find location of closest Re
            [Y, index] = min(diff);
            clear diff
            
            %Interpolate local Cl and Cd based on AoA and Re
            Cl_c(j) = interp1(aoa(:,j),Cl(:,index,j),alpha(j), 'linear', 'extrap');
            Cd_c(j) = interp1(aoa(:,j),Cd(:,index,j),alpha(j), 'linear', 'extrap');
            
            %Determine normal and tangential force coefficients
            C_n(j) = Cl_c(j)*cos(phi(j)) + Cd_c(j)*sin(phi(j));
            C_t(j) = Cl_c(j)*sin(phi(j)) - Cd_c(j)*cos(phi(j));
            
            %Prandlt's Tip Loss Correction
            if sin(phi(j)) < 0.2
               F = 1.0;
            else
               corr = Bl/2*(rotor(end,1) - abe(j))/(abe(j)*sin(phi(j)));
               F = 2*acos(exp(-corr))/pi;
            end
            
            %Glauert Correction for high values of a, Using Spera (1994)
            a_c = 0.2;
             if a(j) > a_c
                 K = 4.*F.*sin(phi(j))^2./(sigma(j)*C_n(j));
                 a_new = (1/2)*(2+K*(1-2*a_c) - sqrt((K*(1-2*a_c)+2)^2 + 4*(K*a_c^2 - 1)));
             else
                a_new = 1./(4*F*sin(phi(j))^2/(sigma(j)*C_n(j))+1);
             end
            
            a_p_new = 1./(4*F*sin(phi(j))*cos(phi(j))/(sigma(j)*C_t(j))-1);
            
            %Relax a-values for stability
            if m < 2
                a_relax = 1;
            elseif m==2
                a_relax = 1;
                a_new = 0.25*a_new + 0.5*a(j) + 0.25*a_prev;
            else
                a_relax = relax;
            end
            
            %Error calculation
            if m == 1
                calc_er = 10;
            else
            calc_er = abs(a_new - a(j));
            end
            if calc_er < tlr
                a(j) = a_new;
                a_p(j) = a_p_new;
                 %disp(['Iteration ' num2str(m)])
                break
            else
                a_prev = a(j);
                a(j) = a_relax*a_new + (1-a_relax)*a(j);
                a_p(j) = a_relax*a_p_new + (1-a_relax)*a_p(j);
            end
            
            clear anew apnew index Y
            m = m + 1;
%             if m == max_iter
%                 disp(['Maximum Iterations Reached ' num2str(a_new - a(j))])
%             end
        end
    
    %Calculate Forces
    P_n(j,iter,i) = 0.5*rho*U_loc(j,iter)^2*abc(j)*C_n(j);
    P_t(j,iter,i) = 0.5*rho*U_loc(j,iter)^2*abc(j)*C_t(j);
     
    l(j) = 0.5*rho*U_loc(j,iter)^2*Cl_c(j)* abe(j);
    d(j) = 0.5*rho*U_loc(j,iter)^2*Cd_c(j)* abe(j);

    p_z(j) = l(j)*cos(phi(j)) + d(j)*sin(phi(j));
    p_z_save(j) = p_z(j); %used to compare values before and after centrifugal force accounted for
    p_y(j) = l(j)*sin(phi(j)) - d(j)*cos(phi(j));
        
 
if iter == 1
   %Calculate induced velocities
    a_star = -W_z(j,iter)/U_loc(j,iter); %Only valid for 0 yaw
    if a_star > 0.5
        a_star = 0.5;
    end
   %Glauert correction
    if a_star <= 0.2
        fg = 1.0;
    else
        fg = 0.2/a_star*(2 - 0.2/a_star);
    end
    
    tau1 = 1.1/(1-1.3*a_star)*abe(end)/U(i); 
    tau2 = (0.39 - 0.26*(abe(j)/abe(end))^2)*tau1;
    W_qsy(j,iter) = -Bl*l(j)*sin(phi(j))/(4*pi*rho*abe(j)*F*((V_rel_z1(j) + fg*W_z(j,iter))^2)^0.5);
    W_qsz(j,iter) = -Bl*l(j)*cos(phi(j))/(4*pi*rho*abe(j)*F*((V_rel_z1(j) + fg*W_z(j,iter))^2)^0.5);
    Hy(j,iter) = W_qsy(j,iter) + 0.6*tau1*(W_qsy(j,iter) - W_qsy(j,iter))/tstep;
    Hz(j,iter) = W_qsz(j,iter) + 0.6*tau1*(W_qsz(j,iter) - W_qsz(j,iter))/tstep;
    W_inty(j,iter) = Hy(j,iter) + (W_inty(j,iter) - Hy(j,iter))*exp(-tstep/tau1);
    W_intz(j,iter) = Hz(j,iter) + (W_intz(j,iter) - Hz(j,iter))*exp(-tstep/tau1);
    W_y(j,iter) = W_inty(j,iter) + (W_y(j,iter) - W_inty(j,iter))*exp(-tstep/tau2);
    W_z(j,iter) = W_intz(j,iter) + (W_z(j,iter) - W_intz(j,iter))*exp(-tstep/tau2); 
elseif iter >=2 && iter <= tmax/tstep
   %Calculate induced velocities
    a_star = -W_z(j,iter-1)/U(i); %Only valid for 0 yaw
    if a_star > 0.5
        a_star = 0.5;
    end
   %Glauert correction
    if a_star <= 0.2
        fg = 1.0;
    else
        fg = 0.2/a_star*(2 - 0.2/a_star);
    end
    
    tau1 = 1.1/(1-1.3*a(j))*abe(end)/U(i);
    tau2 = (0.39 - 0.26*(abe(j)/abe(end))^2)*tau1;
    W_qsy(j,iter) = -Bl*l(j)*sin(phi(j))/(4*pi*rho*abe(j)*F*((V_rel_z1(j) + fg*W_z(j,iter-1))^2)^0.5);
    W_qsz(j,iter) = -Bl*l(j)*cos(phi(j))/(4*pi*rho*abe(j)*F*((V_rel_z1(j) + fg*W_z(j,iter-1))^2)^0.5);
    Hy(j,iter) = W_qsy(j,iter) + 0.6*tau1*(W_qsy(j,iter) - W_qsy(j,iter-1))/tstep;
    Hz(j,iter) = W_qsz(j,iter) + 0.6*tau1*(W_qsz(j,iter) - W_qsz(j,iter-1))/tstep;
    W_inty(j,iter) = Hy(j,iter) + (W_inty(j,iter-1) - Hy(j,iter))*exp(-tstep/tau1);
    W_intz(j,iter) = Hz(j,iter) + (W_intz(j,iter-1) - Hz(j,iter))*exp(-tstep/tau1);
    W_y(j,iter) = W_inty(j,iter) + (W_y(j,iter-1) - W_inty(j,iter))*exp(-tstep/tau2);
    W_z(j,iter) = W_intz(j,iter) + (W_z(j,iter-1) - W_intz(j,iter))*exp(-tstep/tau2);
end

    %Calcualte centrifugal load, N/m
    if j>1
        for jj = 2:numel(radius)
        F_x(jj) = F_x(jj-1) + radius(jj)*speed^2*mass(jj);
        p_cent(jj) = F_x(jj)*(thetaz(jj)-thetaz(jj-1)) - radius(jj)*speed^2*mass(jj)*thetaz(jj);
        end
        p_z(j) = l(j)*cos(phi(j)) + d(j)*sin(phi(j)) + p_cent(j-1);
    end
     
    end

    %Calculate Moment
    for k = 1:(numel(abe)-1)
        A(k,iter,i) = (P_t(k+1,iter,i) - P_t(k,iter,i))/(abe(k+1) - abe(k));
        B(k,iter,i) = (P_t(k,iter,i)*abe(k+1) - P_t(k+1,iter,i)*abe(k))/(abe(k+1) - abe(k));
        M(k,iter,i) = (1/3)*A(k,iter,i)*(abe(k+1)^3 - abe(k)^3) + (1/2)*B(k,iter,i)*(abe(k+1)^2 - abe(k)^2);
    end
    M_tot(iter,i) = Bl*sum(M(:,iter,i));
    
    %Calculate Thrust
    T_n = trapz(abe,P_n(:,iter,i));
    T_tot(iter,i) = Bl*T_n;
    
    %Calculate Coefficient of Power and Thrust
    C_p(iter,i) = (M_tot(iter,i)*speed.*pi./30)./(0.5*rho*U(i)^3*pi*abe(end,1)^2);
    P_out(iter,i) = (M_tot(iter,i)*speed.*pi./30);
    C_thr(iter,i) = T_tot(iter,i)./(0.5*rho*U(i)^2*pi*abe(end)^2);

    t = t+tstep;
    iter = iter+1;
    
end 
disp(['Time completed ' num2str(t)]);
t = 0;
iter = 1;
end

% %Calculate Eigenmodes (first flapwise, first edgewise, second flapwise)
% radius = eim{1};
% twist = eim{3};
% eixx = EIxx;
% eiyy = EIyy;
% [ omegas, uyev, uzev, mbendy, mbendz ] = tara_bladeeigenmode( radius,  twist, eixx, eiyy, mass, pitch, p_y, p_z);
% u1fy = uyev(:,1);
% u1ey = uyev(:,2);
% u2fy = uyev(:,3);
% u3fy = uyev(:,4);
% u1fz = uzev(:,1);
% u1ez = uzev(:,2);
% u2fz = uzev(:,3);
% u3fz = uzev(:,4);
% 
% %Use edgewise mode to calculate blade stiffness
% k_blade(:) = u1ey(:).^2 .* mass(:);
% 
% 
% %Plot Eigenmodes
% hold on
% title('First Flapwise Eigenmode');
% plot(radius,u1fy);
% plot(radius,u1fz);
% legend('u^1^f_y', 'u^1^f_z')
% hold off
% 
% figure(2)
% hold on
% title('First Edgewise Eigenmode');
% plot(radius,u1ey);
% plot(radius,u1ez);
% legend('u^1^e_y', 'u^1^e_z')
% hold off
% 
% figure(3)
% hold on
% title('Second Flapwise Eigenmode');
% plot(radius,u2fy);
% plot(radius,u2fz);
% legend('u^2^f_y', 'u^2^f_z')
% hold off
% 
% figure(4)
% hold on
% title('Third Flapwise Eigenmode');
% plot(radius,u3fy);
% plot(radius,u3fz);
% legend('u^3^f_y', 'u^3^f_z')
% hold off