clear all;
close all;

A =csvread("1541HelixZPos.csv");
%A =csvread("1540Helix.csv");
S = A(:,3);
x0=A(1,1)
y0=A(1,2)
z0=A(1,3)

X= A(:,1)';
Y= A(:,2)';
Z= A(:,3)';


%Initial Condition
    %Initial Guess
      Rperp=22;   %Radius
      Rparl=-10;  %Pitch
    %Produced by Hough Transform (Central Point)
      xp=-54;   
      yp=-110;
    %---<
    theta0=atan2(y0-yp,x0-xp);
    thetaMinus=theta0;
for k=1:length(X)


    theta1(k)=atan2(Y(k)-yp,X(k)-xp); %Imagine this is time elapsed
    theta(k) = theta1(k)-theta0;      %Imagine this is delta-time
    
    if ((theta(k)-thetaMinus)>2)
       %z0=xhatminus(3,k-1)+abs(Rparl*theta(k));     %+20 Because when Theta 2 it Jumped -18
       z0=Z(k)+abs(Rparl*theta(k));
       xhatminus(3,k-1)
       abs(Rparl*theta(k))
    endif
    thetaMinus=theta(k);
    
    
    #Time update
    F=[cos(theta0+theta(k)) 0; sin(theta0+theta(k)) 0; 0 theta(k)];
    xhatminus(:,k)=F*[Rperp Rparl]'+[xp yp z0]';
    
      %MEASUREMENT UPDATE USING Rparl measurement

    %Error Measurement
    errorT(:,k)=[X(k);Y(k);Z(k)]-xhatminus(:,k);
    errorAbs = abs(errorT);
    
endfor
    start=1;
    finnish=length(X);
 
    plot3(xhatminus(3,start:finnish),xhatminus(1,start:finnish),xhatminus(2,start:finnish),'.')
    hold on;
    plot3(Z(start:finnish),X(start:finnish),Y(start:finnish),'r*')
