function [uy,thetay,mbendy,uz,thetaz,mbendz]=deflec2(r,ei1,ei2,twist,pitch,py,pz)
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
for j=ni:-1:2
    tyshear(j-1)=tyshear(j)+0.5*(py(j-1)+py(j))*(r(j)-r(j-1));
    tzshear(j-1)=tzshear(j)+0.5*(pz(j-1)+pz(j))*(r(j)-r(j-1));
end;

mbendy(ni)=0.d0;
mbendz(ni)=0.d0;
for j=ni:-1:2
    mbendy(j-1)=mbendy(j)-tzshear(j)*(r(j)-r(j-1))-(pz(j-1)/6+pz(j)/3)*(r(j)-r(j-1))^2;
    mbendz(j-1)=mbendz(j)+tyshear(j)*(r(j)-r(j-1))+(py(j-1)/6+py(j)/3)*(r(j)-r(j-1))^2;
end;

mprincipal1=1:ni;
mprincipal2=1:ni;
kappa1=1:ni;
kappa2=1:ni;
kappay=1:ni;
kappaz=1:ni;

for j=1:ni;
    mprincipal1(j)=mbendy(j)*cosd(twist(j))-mbendz(j)*sind(twist(j));
    mprincipal2(j)=mbendy(j)*sind(twist(j))+mbendz(j)*cosd(twist(j));
    kappa1(j)=mprincipal1(j)/ei1(j);
    kappa2(j)=mprincipal2(j)/ei2(j);
    kappay(j)=kappa1(j)*cosd(twist(j)+pitch)+kappa2(j)*sind(twist(j)+pitch);
    kappaz(j)=-kappa1(j)*sind(twist(j)+pitch)+kappa2(j)*cosd(twist(j)+pitch);
end;

thetay=1:ni;
thetaz=1:ni;

thetay(1)=0.0;
thetaz(1)=0.0;
for j=1:ni-1;
    thetay(j+1)=thetay(j)+0.5*(kappay(j+1)+kappay(j))*(r(j+1)-r(j));
    thetaz(j+1)=thetaz(j)+0.5*(kappaz(j+1)+kappaz(j))*(r(j+1)-r(j));
end;

uy=1:ni;
uz=1:ni;

uy(1)=0.0;
uz(1)=0.0;
for j=1:ni-1;
    uy(j+1)=uy(j)+thetaz(j)*(r(j+1)-r(j))+(kappaz(j+1)/6.0+kappaz(j)/3.0)*(r(j+1)-r(j))^2;
    uz(j+1)=uz(j)-thetay(j)*(r(j+1)-r(j))-(kappay(j+1)/6.0+kappay(j)/3.0)*(r(j+1)-r(j))^2;
end;
