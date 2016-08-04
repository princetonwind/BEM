function [uy,thetay,mbendy,uz,thetaz,mbendz]=tara_deflec2(r,ei1,ei2,twist,pitch,py,pz, p_y, p_z)
% [uy,thetay,mbendy,uz,thetaz,mbendz]=deflec2(r,ei1,ei2,twist,pitch,py,pz)
% Finds the blade deflection given input loads, by M.O. Hansen %
% See Aerodynamics of Wind Turbines, Ch. 11 for details on eqns%
% r = list of radial stations
% ei1 = E*Ixx; ei2 = E*Iyy;
% twist = local blade twist; pitch = blade pitch angle
% py and pz are the tangential and normal loads, respectively

ni=max(size(r));
% y and z are global coordinates, z is in axial (streamwise) direction
% y is normal to z (looking at x-sec of blade airfoil)
% 1 and 2 are blade-local coordinates or principal axes
tyshear(ni)=0.0;
tzshear(ni)=0.0;

ii = ni;

for ii=ni:-1:2
    tyshear(ii-1)= -p_y(ii)*(r(ii)-r(ii-1)) + tyshear(ii) + 0.5*(py(ii-1)+py(ii))*(r(ii)-r(ii-1));
    tzshear(ii-1)= -p_z(ii)*(r(ii)-r(ii-1)) + tzshear(ii) + 0.5*(pz(ii-1)+pz(ii))*(r(ii)-r(ii-1));
end;

mbendy(ni)=0.d0;
mbendz(ni)=0.d0;
for ii=ni:-1:2
    mbendy(ii-1)=mbendy(ii)-tzshear(ii)*(r(ii)-r(ii-1))-(pz(ii-1)/6+pz(ii)/3)*(r(ii)-r(ii-1))^2;
    mbendz(ii-1)=mbendz(ii)+tyshear(ii)*(r(ii)-r(ii-1))+(py(ii-1)/6+py(ii)/3)*(r(ii)-r(ii-1))^2;
end;

mprincipal1=1:ni;
mprincipal2=1:ni;
kappa1=1:ni;
kappa2=1:ni;
kappay=1:ni;
kappaz=1:ni;

for ii=1:ni;
    mprincipal1(ii)=mbendy(ii)*cosd(twist(ii))-mbendz(ii)*sind(twist(ii));
    mprincipal2(ii)=mbendy(ii)*sind(twist(ii))+mbendz(ii)*cosd(twist(ii));
    kappa1(ii)=mprincipal1(ii)/ei1(ii);
    kappa2(ii)=mprincipal2(ii)/ei2(ii);
    kappay(ii)=kappa1(ii)*cosd(twist(ii)+pitch)+kappa2(ii)*sind(twist(ii)+pitch);
    kappaz(ii)=-kappa1(ii)*sind(twist(ii)+pitch)+kappa2(ii)*cosd(twist(ii)+pitch);
end;

thetay=1:ni;
thetaz=1:ni;

thetay(1)=0.0;
thetaz(1)=0.0;
for ii=1:ni-1;
    thetay(ii+1)=thetay(ii)+0.5*(kappay(ii+1)+kappay(ii))*(r(ii+1)-r(ii));
    thetaz(ii+1)=thetaz(ii)+0.5*(kappaz(ii+1)+kappaz(ii))*(r(ii+1)-r(ii));
end;

uy=1:ni;
uz=1:ni;

uy(1)=0.0;
uz(1)=0.0;
for ii=1:ni-1;
    uy(ii+1)=uy(ii)+thetaz(ii)*(r(ii+1)-r(ii))+(kappaz(ii+1)/6.0+kappaz(ii)/3.0)*(r(ii+1)-r(ii))^2;
    uz(ii+1)=uz(ii)-thetay(ii)*(r(ii+1)-r(ii))-(kappay(ii+1)/6.0+kappay(ii)/3.0)*(r(ii+1)-r(ii))^2;
end;
