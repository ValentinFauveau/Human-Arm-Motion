%Simulation of Human Arm Motion

clc;
close all;
clear all;

%Common values
%References: 
%_Rudolf s Drillis, Renato Contini and Maurice Bluestein, 1966,"Body Segment Parameters"
%_John T. McConville, Thomas D. Churchili, 1981, "Anthropometric relationships of body and body segment moments of inertia"

%Mean Height is 177.5 cm

%% Inputs
%Initial position
qs=-pi/6;
qe=pi/4;
qw=pi/8;

%Final position
xf=-20; %In cm
yf=40; %In cm

%% 

%Length of each limb
la=36.4; %In cm
lf=29.9; %In cm
lh=20.3; %In cm

%Mass of each limb
Ma=2.07; %In g
Mf=1.16; %In g
Mh=0.54; %In g

%The center of mass position of each limb
%Distance from proximal joint segment length =1
lma=0.427; %In %
lmf=0.417; %In %
lmh=0.361; %In %

%Position of xi and yi for each joint 
%Shoulder
xs=0;
ys=0;

%Elbow
xe=la*cos(qs);
ye=la*sin(qs);

%Wrist
xw=xe+lf*cos(qs+qe);
yw=ye+lf*sin(qs+qe);

%Finger tips
x=xw+lh*cos(qs+qe+qw);
y=yw+lh*sin(qs+qe+qw);

%Inverse kinematics
teta=qs+qe+qw;

xw2=x-lh*cos(teta);
yw2=y-lh*sin(teta);

%Now we find the angle alfa and l
% l=sqrt(xw^2+yw^2);
l=sqrt((x-lh*cos(teta))^2+(y-lh*sin(teta))^2);
alfa=acos((x-lh*cos(teta))/l);%Angle between horizontal and l

%Cosine law
gama=acos((lf^2-la^2-l^2)/(-2*la*l));%Angle between arm and l

beta=alfa-gama;

xe2=la*cos(beta);
ye2=la*sin(beta);

xs2=0;
ys2=0;

%Jacobian
%x=la*cos(qs)+lf*cos(qs+qe)+lh*cos(qs+qe+qw)
%y=la*sin(qs)+lf*sin(qs+qe)+lh*sin(qs+qe+qw)
%teta=qs+qe+qw

dxdqs=la*(-sin(qs))+lf*(-sin(qs+qe))+lh*(-sin(qs+qe+qw));
dxdqe=lf*(-sin(qs+qe))+lh*(-sin(qs+qe+qw));
dxdqw=lh*(-sin(qs+qe+qw));

dydqs=la*cos(qs)+lf*cos(qs+qe)+lh*cos(qs+qe+qw);
dydqe=lf*cos(qs+qe)+lh*cos(qs+qe+qw);
dydqw=lh*cos(qs+qe+qw);

dtdqs=1;
dtdqe=1;
dtdqw=1;

J=[dxdqs dxdqe dxdqw; dydqs dydqe dydqw; dtdqs dtdqe dtdqw];

%% 
xi=x;
yi=y;

t=0:0.01:5;
d=5;

%Function describing the minimum jerk movement
xt=zeros(2,length(t));
xt(1,:)=xi+(xf-xi)*(6*(t/d).^5-15*(t/d).^4+10*(t/d).^3);
xt(2,:)=yi+(yf-yi)*(6*(t/d).^5-15*(t/d).^4+10*(t/d).^3);

%Different plots
figure;
plot(t,xt(1,:));
title('Evolution of the hand in X direction');
xlabel('Time (s)');
ylabel('X position');

figure;
plot(t,xt(2,:));
title('Evolution of the hand in Y direction');
xlabel('Time (s)');
ylabel('Y position');

figure;
plot(xt(1,:),xt(2,:));
title('Y evolution vs X evolution');
xlabel('X position');
ylabel('Y position');

dxt=diff(xt(1,:));
t=t(1:length(t)-1);
figure;
plot(t,dxt);
title('Angular velocity');
xlabel('Time (s)');
ylabel('Angular velocity (rad/s)');

%Calculation of the joint angle and angular velocity that minimize the
%angular velocity
%We know that dq=J+*dx
%We also know the variation of the x points
%We also know the Jacobian at the first instance

dyt=diff(xt(2,:));

dxdt=zeros(2,length(t));
dxdt=[dxt;dyt];

Q=zeros(3,length(t));
J=J([1 2],:);

%This loop is used to calculate the qs, qe and qw at each instance
%Algorithm: 1)Calculate the Jcross
%           2)Calculate the Qdot
%           3)Calculate the new qs, qe and qw
%           4)Calculate the new J
%           5)Repeat all the steps to all the moments of t.
for i=1:length(t)
    
    if i~=1
        
        dxdqs=la*(-sin(qs))+lf*(-sin(qs+qe))+lh*(-sin(qs+qe+qw));
        dxdqe=lf*(-sin(qs+qe))+lh*(-sin(qs+qe+qw));
        dxdqw=lh*(-sin(qs+qe+qw));
        
        dydqs=la*cos(qs)+lf*cos(qs+qe)+lh*cos(qs+qe+qw);
        dydqe=lf*cos(qs+qe)+lh*cos(qs+qe+qw);
        dydqw=lh*cos(qs+qe+qw);
        
        J=[dxdqs dxdqe dxdqw; dydqs dydqe dydqw];
    end
    Jcross=J'*inv(J*J');
    Q(:,i)=Jcross*dxdt(:,i);
    Q(1,i)=qs+Q(1,i);
    Q(2,i)=qe+Q(2,i);
    Q(3,i)=qw+Q(3,i);
    qs=Q(1,i);
    qe=Q(2,i);
    qw=Q(3,i);
end

%% 

%To generate the movie we first calculate the coordinates for each joint at
%each instance.
Xs=zeros(1,length(t));
Ys=zeros(1,length(t));

Xe=zeros(1,length(t));
Ye=zeros(1,length(t));

Xw=zeros(1,length(t));
Yw=zeros(1,length(t));

Xe=la*cos(Q(1,:));
Ye=la*sin(Q(1,:));

Xw=Xe+lf*cos(Q(1,:)+Q(2,:));
Yw=Ye+lf*sin(Q(1,:)+Q(2,:));

X=Xw+lh*cos(Q(1,:)+Q(2,:)+Q(3,:));
Y=Yw+lh*sin(Q(1,:)+Q(2,:)+Q(3,:));

%Now we generate the movie
figure;
mx=15*cos(linspace(0,2*pi,500))-20;
my=20*sin(linspace(0,2*pi,500))+20;
for i=1:length(t)
    clf;
    plot([0 Xe(i) Xw(i) X(i)] , [0 Ye(i) Yw(i) Y(i)]); %Plot the moving arm
    hold on;
    plot([0 -40],[0 0]); %Plot shoulders
    plot([0 0],[0 -80]); %Plot the body
    plot([-40 -40],[0 -80]); %Plot the body
    plot([0 -40],[-80 -80]); %Plot the body
    plot([-40 -50],[0 -80]); %Plot the other arm
    plot(mx,my); %Plot the head
    axis([-(la+lf+lh)  la+lf+lh -(la+lf+lh)  la+lf+lh]);
    title('Movie of the moving arm');
    hold on;
    pause(0.01);
end

