%Model V27 BEM Code
%Tara Nealon

clc, clear;

%Airfoil Input Data: Column 1 = AoA, Row 1 = Re, Rest of Columns = Cd or CL
aerodata_CD(:,:,1) = importdata('C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\Circle_CD.txt', '\t');
aerodata_CL(:,:,1) = importdata('C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\Circle_CL.txt', '\t');
aerodata_CD(:,:,2) = importdata('C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\N63214_CD.txt', '\t');
aerodata_CL(:,:,2) = importdata('C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\N63214_CL.txt', '\t');
aerodata_CD(:,:,3) = importdata('C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\N63218_CD.txt', '\t');
aerodata_CL(:,:,3) = importdata('C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\N63218_CL.txt', '\t');
aerodata_CD(:,:,4) = importdata('C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\N63222_CD.txt', '\t');
aerodata_CL(:,:,4) = importdata('C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\N63222_CL.txt', '\t');
aerodata_CD(:,:,5) = importdata('C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\N63228_CD.txt', '\t');
aerodata_CL(:,:,5) = importdata('C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\N63228_CL.txt', '\t');

%data_CD(:,1,:) = aerodata_CD(2:62,1,:);

%Choose a Reynolds Number for Data
Re = 1500000;
[k,j] = find(aerodata_CD(:,:,1) == Re);
data_CD(:,2,1) = aerodata_CD(2:62,j,1);
data_CD(:,2,2) = aerodata_CD(2:62,j,2);
data_CD(:,2,3) = aerodata_CD(2:62,j,3);
data_CD(:,2,4) = aerodata_CD(2:62,j,4);
data_CD(:,2,5) = aerodata_CD(2:62,j,5);
data_CL(:,2,1) = aerodata_CL(2:62,j,1);
data_CL(:,2,2) = aerodata_CL(2:62,j,2);
data_CL(:,2,3) = aerodata_CL(2:62,j,3);
data_CL(:,2,4) = aerodata_CL(2:62,j,4);
data_CL(:,2,5) = aerodata_CL(2:62,j,5);
data_CD(:,1,:) = aerodata_CD(2:62,1,:);
data_CL(:,1,:) = aerodata_CL(2:62,1,:);

%Define Geometry (9 blade section): radius (m), chord (m), twist (deg), airfoil type (1 = %circle, 2= N63214, 3 = N63218, 4 = N63222, 5 = N63228)
geometry = [0.0015 0.00296 13 1 ; 0.0074 0.0030 13 1; 0.0148 0.0096 13 5; 0.0296 0.0086 10.7 5; 0.0444 0.0076 8.5 5; 0.0593 0.0065 6.2 4; 0.0741 0.0055 4 3; 0.0889 0.0044 1.7 2; 0.1000 0.0037 0 2];

%Definte Constants
w_1 = 43*2*pi/60; %Rotational speed in rad/s, generator 1, factor of 135 for to convert to model RPM
w_2 = 33*2*pi/60; %Rotational speed in rad/s, generator 2, factor of 135 for to convert to model RPM
B = 3; %Number of blades
V_cut_in = 3.5; %Cut-in Wind Speed, m/s
V_cut_out = 25; %Cut-off Wind Speed, m/s
V_rated = 14; %Rated Wind Speed, m/s
P_rated = 225000; %Rated Power, kW
rho = 1.225%*210; %kg/m3, multiplied by 210 to convert to model density
mu = (1.79*10^-5)%/1.2; %kg/m s, divided by 1.2 to convert to model viscosity
a_c = 0.2; %Critical axial induction factor 
relax = 0.5;
tolerance = 1E-6;

%Initialize Parameters
a = zeros(9,1);
a_new = zeros(9,1);
a_prime = zeros(9,1);
a_prime_new = zeros(9,1);
wind_type = 1; %1 = cst
theta_p = 0; %Pitch angle, degrees
phi = zeros(9,1);
theta = zeros(9,1);
alpha = zeros(9,1);
V_rel = zeros(9,2);
W = zeros(9,2);
Cl = zeros(9,1);
Cd = zeros(9,1);
l = zeros(9,1);
d = zeros(9,1);
V_rel_mag = zeros(9,1);
p = zeros(9,2);
W_qs = zeros(9,2);
W = zeros(9,2);
W_old = zeros(9,2);
W_qs_old = zeros(9,2);
relax = 0.1;
P=zeros(22,1);
P_r=zeros(22,1);
x=1;

%Define Wind Inlet
if wind_type == 1 %constant inlet velocity
    V_o(1) = 0; %y-direction, m/s
    V_o(2) = 3.5; %z-direction, m/s
    V_o_mag = sqrt(V_o(1)^2 + V_o(2)^2);
end

error = 10;
iter = 0;
%Begin BEM Simulation

while V_o(2) <= 25
    while error > 0.0001
        for i = 1:9
            %Compute Relative Velocity Components
            V_rel(i,1) = V_o(1) - w_1*geometry(i,1) + W_old(i,1); %y-direction, m/s
            V_rel(i,2) = V_o(2) + W_old(i,2); %z-direction, m/s
            
            %Compute Local Flow Angle phi
            phi(i) = atan((1-a(i)*V_o(2))/(1+a_prime(i)*w_1*geometry(i,1)));
            
            %FIX PHI
            
            %Compute Local AoA
            alpha(i) = 180/pi*phi(i) - geometry(i,3) - theta_p; %AoA, degrees
            
            %Look Up Cd and Cl
            
            if i <=2
                type = 1;
                Cd(i) = interp1(data_CD(:,1,type),data_CD(:,2,type),alpha(i),'linear','extrap');
                Cl(i) = interp1(data_CL(:,1,type),data_CL(:,2,type),alpha(i),'linear','extrap');
            elseif i>=3 && i<= 5
                type = 5;
                Cd(i) = interp1(data_CD(:,1,type),data_CD(:,2,type),alpha(i),'linear','extrap');
                Cl(i) = interp1(data_CL(:,1,type),data_CL(:,2,type),alpha(i),'linear','extrap');
            elseif i==6
                type = 4;
                Cd(i) = interp1(data_CD(:,1,type),data_CD(:,2,type),alpha(i),'linear','extrap');
                Cl(i) = interp1(data_CL(:,1,type),data_CL(:,2,type),alpha(i),'linear','extrap');
            elseif i == 7
                type = 3;
                Cd(i) = interp1(data_CD(:,1,type),data_CD(:,2,type),alpha(i),'linear','extrap');
                Cl(i) = interp1(data_CL(:,1,type),data_CL(:,2,type),alpha(i),'linear','extrap');
            elseif i>=8 && i<= 9
                type = 2;
                Cd(i) = interp1(data_CD(:,1,type),data_CD(:,2,type),alpha(i),'linear','extrap');
                Cl(i) = interp1(data_CL(:,1,type),data_CL(:,2,type),alpha(i),'linear','extrap');
            end
            
            %Calculate relative velocity
            V_rel_mag(i) = sqrt(V_rel(i,1)^2 + V_rel(i,2)^2);
            
            %Calculate Lift and Drag Forces per Length
            l(i) = 0.5*rho*V_rel_mag(i)^2 * Cl(i) * geometry(i,2);
            d(i) = 0.5*rho*V_rel_mag(i)^2 * Cd(i) * geometry(i,2);
            
            %Tangential and Normal Forces per Length
            p(i,1) = l(i)*sin(phi(i)) - d(i)*cos(phi(i)); %tangential force (y-direction)
            p(i,2) = l(i)*cos(phi(i)) - d(i)*sin(phi(i)); %normal force (z-direction)
            
            if i == 9
                p(i,1) = 0;
                p(i,2) = 0;
            end

            %Normal and Tangential Coefficients
            C_n(i) = p(i,2)/(1/2*rho*V_rel_mag(i)^2*geometry(i,2));
            C_t(i) = p(i,1)/(1/2*rho*V_rel_mag(i)^2*geometry(i,2));
            
            %Solidity
            sigma(i) = geometry(i,2)*B/(2*pi*geometry(i,1));
            
            %a_prime(i) = 1/(4*sin(phi(i))*(cos(phi(i)))/ (sigma(i) * C_t(i)) - 1);
            
            %Valid for zero yaw
            a_star = -W_old(i,2)/sqrt(V_o(1)^2 + V_o(2)^2);
            
%             %Glauert Correction, a_c = critical axial induction factor
%             if a_star<a_c
%                 fg = 1.0;
%             else
%                 fg = a_c/a_star*(2 - a_c/a_star);
%             end
%             
            %Compute term for denominator of induced wind calculation
            %glauert_coeff = sqrt(V_o(1)^2 + (V_o(2) + fg*W_old(i,2))^2);
            
            %Tip Loss Correction Factor
            if sin(phi(i)) < 0.02
                F = 1;
            else
                F = 2/pi*acos(exp(-B*(geometry(9,1)-geometry(i,1)+0.001)/(2*geometry(i,1)*sin(phi(i)))));
            end
            
            %Glauert Correction, Spera (1994)
            if a(i) > 0.2
                corr = 4.*F.*sin(phi(i))^2./(sigma(i)*C_n(i));
                a_new(i) = 0.5*(2 + corr*(1-2*0.2) - sqrt((corr*(1-2*0.2) + 2)^2) + 4*(corr*0.2^2 - 1));
            else
                a_new(i) = 1./((4*F*sin(phi(i)).^2)/(sigma(i)*C_n(i))+1);
            end
            
            a_prime_new(i) = 1./((4*F*sin(phi(i))*cos(phi(i)))/(sigma(i)*C_t(i))-1);
                            
            if iter < 2
                a_relax(i) = 1;
            elseif iter == 2
                a_relax(i) = 1;
                a_new(i) = 0.25*a_new(i) + 0.5*a(i) + 0.25*a_prev(i);
            else
                a_relax(i) = relax;
            end
                        
            if abs(a_new - a(i)) < tolerance
                %a_prev(i) = a(i);
                a(i) = a_new(i);
                a_prime(i) = a_prime_new;
            else
                a_prev(i) = a(i);
                a(i) = a_relax(i)*a_new(i) + (1-a_relax(i))*a(i);
                a_prime(i) = a_relax(i)*a_prime_new(i) + (1-a_relax(i))*a_prime(i);
            end
            
            
            %Begin Dynamic Wake
%             W_new(i,1) = (-B*l(i)*sin(phi(i)))/(4*pi*rho*geometry(i,1)*F*glauert_coeff);
%             W_new(i,2) = (-B*l(i)*cos(phi(i)))/(4*pi*rho*geometry(i,1)*F*glauert_coeff);
%             
%             error1 = max(abs(W_new(i,1)-W_old(i,1)));
%             error2 = max(abs(W_new(i,2)-W_old(i,2)));
%             error = max(error1,error2);
%             
%             Reset old values to a relaxed new value
%             W_old(i,1) = (1-relax)*W_old(i,1) + relax*W_new(i,1);
%             W_old(i,2) = (1-relax)*W_old(i,2)+ relax*W_new(i,2);
%             
            %Begin Moment and Thrust Calculations
            if i < 9
                A(i) = (p(i+1,1) - p(i,1))/(geometry(i+1,1)-geometry(i,1));
                A_tan(i) = (p(i+1,2) - p(i,2))/(geometry(i+1,1)-geometry(i,1));
                B2(i) = (p(i,1)*geometry(i+1,1) - p(i+1,1)*geometry(i,1))/(geometry(i+1,1) - geometry(i,1));
                B2_tan(i) = (p(i,2)*geometry(i+1,1) - p(i+1,2)*geometry(i,1))/(geometry(i+1,1) - geometry(i,1));
                M(i) = (1/3)*A(i)*(geometry(i+1,1)^3 - geometry(i,1)^3) + (1/2)*B2(i)*(geometry(i+1,1)^2 - geometry(i,1)^2);
                T(i) = 1/2*A_tan(i)*(geometry(i+1,1)^2 - geometry(i,1)^2) + B2_tan(i)*(geometry(i+1,1) - geometry(i,1));
            end 
        end
        
        %Moment, Power, and Thrust Calculations
        M_tot(x) = sum(M)*B;
        P(x) = M_tot(x)*w_1/1000; %Calculate Power in kW
        P_r(x) = 2.25;  %Rated Power in kW
        T_tot(x) = sum(T)*B/1000;  %Thrust
        C_t(x) = T_tot(x)*1000/(1/2*rho*pi*geometry(9,1)^2*V_o(2)^2); %Coefficient of Thrust
        C_p(x) = (M_tot(x)*w_1*1000)./(0.5*rho*V_o_mag^3*geometry(9,1)^2);
        %C_p(x) = P(x)*1000/(1/2*rho*pi*geometry(9,1)^2*V_o(2)^3); %Coefficient of Power
        
        iter = iter + 1;
        error
    end
    V_o(2) = V_o(2)+1
    x = x+1;
    error = 10;
end

hold on
plot(5:26,P)
plot(5:26,P_r)



