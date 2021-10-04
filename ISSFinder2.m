    function[Position1,Velocity1,i,omega,e,p,a,w,v,v2,Position2,Velocity2,IterationHist]=ISSFinder2(Range,Range_Rate,Latitude,Longitude,Theta_g0,Year,Month1,Day1,Month2,Day2,Hour2,Minute2,TOF)
t1 = datetime(Year,Month1,Day1);
t1 = day(t1,'dayofyear');
t2 = datetime(Year,Month2,Day2);
t2 = day(t2,'dayofyear');
dt=t2-t1;
fraction = ((Hour2.*60)+Minute2)./1440;
daycount = dt + fraction;
Theta_g1 = Theta_g0 + (1.00273779093.*360.*daycount);
Theta_g2 = Theta_g1./360;
Theta_g = Theta_g1 - (360.*(floor(Theta_g2)));
Theta = Theta_g - Longitude;
rh = [Range(1), Range(2), Range(3) + 1];
rr = Range_Rate;
D11= sind(Latitude).*cosd(Theta);
D12= -1.*sind(Theta);
D13= cosd(Latitude).*cosd(Theta);
D21= sind(Latitude).*sind(Theta);
D22= cosd(Theta);
D23= cosd(Latitude).*sind(Theta);
D31= -1.*cosd(Latitude);
D32= 0;
D33= sind(Latitude);
X= D11.*rh(1) + D12.*rh(2) + D13.*rh(3);
Y= D21.*rh(1) + D22.*rh(2) + D23.*rh(3);
Z= D31.*rh(1) + D32.*rh(2) + D33.*rh(3);
Position1 = [X,Y,Z];
Xv= D11.*rr(1) + D12.*rr(2) + D13.*rr(3);
Yv = D21.*rr(1) + D22.*rr(2) + D23.*rr(3);
Zv = D31.*rr(1) + D32.*rr(2) + D33.*rr(3);
range_rate = [Xv,Yv,Zv];
wcrossr = [-1.*.05883.*Y,.05883.*X,0];
Velocity1 = range_rate + wcrossr;
P=Position1;
V=Velocity1;
h_vector=[(P(2).*V(3))-(P(3).*V(2)),(P(3).*V(1))-(P(1).*V(3)),(P(1).*V(2))-(P(2).*V(1))];
h_mag=sqrt((h_vector(1).^2+h_vector(2).^2+h_vector(3).^2));
i=acosd(h_vector(3)./h_mag);
n_vector=[-1.*h_vector(2),h_vector(1),0];
n_mag=sqrt((n_vector(1).^2+n_vector(2).^2));
omega=acosd(n_vector(1)./n_mag);
if omega < 180 & n_vector(2) > 0
    omega = omega;
elseif omega > 180 & n_vector(2) > 0
    omega = 360 - omega;
elseif omega < 180 & n_vector(2) < 0
    omega = 360 - omega;
elseif omega > 180 & n_vector(2) < 0
    omega = omega;
end
r_mag=sqrt((P(1).^2+P(2).^2+P(3).^2));
v_mag=sqrt((V(1).^2+ V(2).^2+ V(3).^2));
rdotv=(P(1).*V(1))+(P(2).*V(2))+(P(3).*V(3));
eterm1=v_mag.^2-(1./r_mag);
evec1=[eterm1.*P(1),eterm1.*P(2),eterm1.*P(3)];
evec2=[rdotv.*V(1),rdotv.*V(2),rdotv.*V(3)];
evec=evec1-evec2;
e=sqrt((evec(1).^2+evec(2).^2+evec(3).^2));
p=h_mag.^2;
ndote=(n_vector(1).*evec(1))+(n_vector(2).*evec(2))+(n_vector(3).*evec(3));
w=acosd(ndote./(n_mag.*e));
if w < 180 & evec(3) > 0
    w = w;
elseif w > 180 & evec(3) > 0
    w = 360 - w;
elseif w < 180 & evec(3) < 0
    w = 360 - w;
elseif w > 180 & evec(3) < 0
    w = w;
end
edotr=(evec(1).*P(1))+(evec(2).*P(2))+(evec(3).*P(3));
v=acosd(edotr./(e.*r_mag));
if v < 180 & rdotv > 0
    v = v;
elseif v > 180 & rdotv > 0
    v = 360 - v;
elseif v < 180 & rdotv < 0
    v = 360 - v;
elseif v > 180 & rdotv < 0
    v = v;
end
val1=tand(v/2);
val2=sqrt((1+e)/(1-e));
M0=atand(val1/val2);
M0=M0.*(pi/180);
a=p/(1-e.^2);
a2=a.*6378.145;
Period = (6.28).*(sqrt((a.^3)/398600));
peripasses = TOF./Period;
peripasses = floor(peripasses);
M = ((sqrt(398600./(a2.^3))).*TOF)-(6.28.*peripasses) + M0;
E=.5;
dM=1;
num=0;
b=[];
c=[];
d=[];
f=[];
while abs(dM) > 1.e-7
    b=[b;E];
    dM = M - (E-(e.*sin(E)));
    c=[c;dM];
    dMdE = 1-(e.*cos(E));   
    d=[d;dMdE];
    E = E + (dM./dMdE);
    f=[f;E];
    num=num+1;
end
iteration=[0:(num-1)];
iteration=transpose(iteration);
IterationHist=table(iteration,b,c,d,f);
IterationHist.Properties.VariableNames = {'i', 'Ei','dM','dMdE','EiNext'};
v2=acos((cos(E)-e)./1-(e.*cos(E)));
v2=v2.*(180/pi);
rw=[r_mag.*cosd(v2),r_mag.*sind(v2)];
vw=[sqrt(1/p).*-1.*sind(v2),sqrt(1/p).*(e+cosd(v2))];
R11=(cosd(omega).*cosd(w))-(sind(omega).*sind(w).*cosd(i));
R12=(-1.*cosd(omega).*sind(w))-(sind(omega).*cosd(w).*cosd(i));
R13=sind(omega).*sind(i);
R21=(sind(omega).*cosd(w))+(cosd(omega).*sind(w).*cosd(i));
R22=(-1.*sind(omega).*sind(w))+(cosd(omega).*cosd(w).*cosd(i));
R23=-1.*cosd(omega).*sind(i);
R31=sind(w).*sind(i);
R32=cosd(w).*sind(i);
R33=cosd(i);
Rx=(rw(1).*R11)+(rw(2).*R12);
Ry=(rw(1).*R21)+(rw(2).*R22);
Rz=(rw(1).*R31)+(rw(2).*R32);
Vx=(vw(1).*R11)+(vw(2).*R12);
Vy=(vw(1).*R21)+(vw(2).*R22);
Vz=(vw(1).*R31)+(vw(2).*R32);
Position2 = [Rx Ry Rz];
Velocity2 = [Vx Vy Vz];
end