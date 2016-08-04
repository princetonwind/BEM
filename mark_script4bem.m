% QBLADE COMPARISONS %
clear all
close all
clc
    T = 20;
    P = 103100;
    speed = 43;
    pitch = 0;
    U = 3.5:0.5:25;
    nb = 3;
    NE = 9;
    relax = 0.5;
    scale = 1.0;
%%    
    %Rotor geometry file location
    %change this to location of rotor file
    rtype = 'C:\Users\tnealon\Documents\Princeton Research\V27 Data\V27_Full-Scale_BEM_Geometry.txt';
    %Format input file
    fid = fopen(rtype,'r');
    rotorgem = textscan(fid,'%f %f %f %s','Delimiter',' ','MultipleDelimsAsOne',1,'Headerlines',1);
    fclose(fid);
    rotor = [rotorgem{1} rotorgem{2} rotorgem{3}];
    foil = rotorgem{4};
    
%%  

    for ii = 1:numel(U)
    disp(['Currrent Velocity ' num2str(U(ii))])


   [ rho, mu, summ(ii,:), bldout(:,:,ii) ] = mark_bem360(rotor, foil, T, P, speed, pitch , U(ii) , nb, NE) ;

   
    end
    
axisfsize = 30;
    %Plotting power Curves%
        figure(12)
        [C,h] = contourf(U,pitch,summ(:,8)');
        set(h,'LevelList',[0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5])
        grid on
        clabel(C,h,'manual','Fontsize',axisfsize,'Color','w','FontWeight','bold')
        ylabel('Blade Pitch (degrees)','FontUnits','points','FontSize',axisfsize)
        xlabel('Free-stream Velocity','FontUnits','points','FontSize',axisfsize)
        set(gca,'FontName','Times','FontUnits','points','FontSize',28)
        
    %Interpolating to find pitch angle%
%     pitches = zeros(numel(U),1);
%     for m = 1:numel(U)
%         if U(m) >= 14
%         pitches(m) = interp1(summ(:,8),pitch,2);
%         end
%     Ctf(m) = interp1(pitch,Ctout(m,:),pitches(m));
%     Cpf(m) = interp1(pitch,Cpout(m,:),pitches(m));
%     end
%     
    %figure(13)
%    plot(U,pitches,'r','Linewidth',2.5)
   % xlabel('Freestream Velcity (m/s)','FontUnits','points','FontSize',axisfsize)
    %ylabel('Pitch (deg)','FontUnits','points','FontSize',axisfsize)
   % set(gca,'FontName','Times','FontUnits','points','FontSize',28)
    
   % figure(14)
    %[AX,H1,H2] = plotyy(U,Cpf,U,Ctf);
   % xlabel('Freestream Velcity (m/s)','FontUnits','points','FontSize',axisfsize)
%    ylabel(AX(1),'C_p','FontUnits','points','FontSize',axisfsize)

    %ylabel(AX(2),'C_t','FontUnits','points','FontSize',axisfsize)
    %set(AX(1),'FontName','Times','FontUnits','points','FontSize',28)
   %set(AX(2),'FontName','Times','FontUnits','points','FontSize',28)
   % set(H1,'Linewidth',2)
    %set(H2,'Linewidth',2)
    
