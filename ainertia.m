function [ Ixx, Iyy, Atot,ymin,ymax] = ainertia( foil, c)
%This function finds the area moment of inertia for airfoils
% foil = name of airfoil
% c = chord of foil
close all
npts = 500;    %Number of points to interpolate grid
csys = [0.25 0].*c; %Loc. of inertial coordinate system, default at 1/4 c

% disp('Airfoil format: x = [1 to 0 to 1]')
% disp(' and corresponding y points, Ixx & Iyy not correct if otherwise!')

rfolder = 'C:\Users\tnealon\Documents\Princeton Research\V27 Data\Airfoil Data\';
% Load in airfoil Data %
data = dlmread([rfolder char(foil) '.dat']);
x = data(:,1).*c; y = data(:,2).*c;
clear data

hal = find(x == 0, 1);
if isempty(hal) == 1
    warning('Format of airfoil input is not correct! Need x = 0.0 point!')
elseif hal ~= (numel(x) + 1)/2 ;
    warning(['Airfoil must decrease in x to zero then increase'...
        '(top then bottom of airfoil in y)'])
end

x = x - csys(1); y = y - csys(2); %Shift to desired chord location

Atot = 0;
% Interpolate to grid for more accurate answer %
xu = linspace(x(1),x(hal),npts);
yu = interp1(x(1:hal),y(1:hal),xu,'pchip');

xl = linspace(x(hal),x(end),npts);
yl = interp1(x(hal:end),y(hal:end),xl,'pchip');
    plot(x,y,'r-',xu,yu,'bo',xl,yl,'gx')
    hold on
% Approximate Area of inertia %
    for i = 2:numel(xu)
        %Find rectangular section for top%
        b = xu(i-1) - xu(i);
        hu = (yu(i) + yu(i-1)) / 2; % Use mean value for height
        hl = (yl(i) + yl(i-1)) / 2;
        h = hu - hl;                %subtract upper from lower bound
        Ixr(i) = b * h^3 / 12;
        Iyr(i) = b^3 * h / 12;
        %Parallel Axis theorem to shift element to correct location%
        dy = (hu + hl) / 2; dx = (xu(i) + xu(i-1)) / 2;
        A = b * h;
         Ixr(i) = Ixr(i) + A * dy ^ 2;
         Iyr(i) = Iyr(i) + A * dx ^ 2;
        Atot = Atot + A;  
        plot([xu(i-1) xu(i)],[hu hu],'r*',[xl(i-1) xl(i)],[hl hl],'ms',...
            [dx dx],[yu(i) yu(i-1)],'kx')
    end
    
    Ixx = sum(Ixr);
    Iyy = sum(Iyr);
    %Shift total second moment of area to desired location%

    
end

