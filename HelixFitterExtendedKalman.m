function [Ak,P,Pm,residue] = HelixFitterExtendedKalman(M =  csvread("1541HelixZPos.csv"),Pc =[-55 -110]', r= 22, tanLamda = (-60 / (2 * pi) )/r,processError=[5 5 5 5 5]',measError=[0.1 0.1 0.1]')

% INPUT
%
% M  : measurement points should be Z ordered
% Pc : central point
% r  : radius of helix (guess)
% tanLamda : tanLamda of helix (guess) Rp = R*tanLamda
% processErrorVariant : initial process error Variant
% measErrorVariant : initial measurement error Variant
%
% Output
% Ak : process model for each pivot A(:,i) = [dr_i theta_i dz_i R tanLamda]
% P  : measurement points
% Pm : estimation points
% residue : distance between P and Pm

xc = Pc(1);
yc = Pc(2);

% first guess
r = r;
tanLamdaFirst = tanLamda;

finish = 75;


% M  = csvread("1541HelixZPos.csv");



x0 =M(1,1);
y0 =M(1,2);
z0 =M(1,3);

X= M(:,1)';
Y= M(:,2)';
Z= M(:,3)';


theta0 = atan2(y0 -yc  ,x0 - xc);


if (theta0 > pi)
	theta0 = theta0 - 2*pi;
endif
if (theta0 < -pi)
	theta0 = 2*pi + theta0;
endif
r0(1) = sqrt((y0 - yc)*(y0 - yc) + (x0 - xc)*(x0 - xc));
dR = sqrt((y0 - yc)*(y0 - yc) + (x0 - xc)*(x0 - xc)) -   r;  
z0 =M(1,3);

Theta(1) = theta0;
Dr(1) = dR;
%eA = [0.1 0.01 0 0 0]';
%Ck = cov(eA * eA');
Ck = zeros(5,5*length(X));
Fk = zeros(5,5*length(X));
%Ck
Ck(1:5,1:5) = processError *  processError' .* eye(5);
Xm(1) = x0;
Ym(1) = y0;
Zm(1) = z0;


Vk = measError *  measError' .* eye(3);
Gk = Vk^-1;
Rk(1) = r;
tanLamdak(1) = tanLamda;
Ak(:,1) = [Dr(1) theta0 0 r tanLamda]';

error(1) = 0;

for i=2:length(X)
 Ck1 = Ck(1:5,(i-2)*5+1:(i-1)*5) ;
 Ak1 = Ak(:,i-1);
 
 
 Theta(i) = atan2(M(i,2) - yc,M(i,1)-xc);
 % normalize Theta  
 if (Theta(i) > pi)
	Theta(i) = Theta(i) - 2*pi;
 endif
 if (Theta(i) < -pi)
	Theta(i) = 2*pi + Theta(i);
 endif

 
 x =M(i,1);
 y =M(i,2);
 z =M(i,3);
 
 %% new 
 theta = Theta(i);


 deltaTheta = theta - theta0;
 % normalize deltaTheta  
 if (deltaTheta > pi)
	deltaTheta = deltaTheta- 2*pi;
 endif
 if (deltaTheta < -pi)
	deltaTheta = 2*pi + deltaTheta;
 endif


 R0(i) = sqrt((y - yc)*(y - yc) + (x - xc)*(x - xc));
 Dr(i) = sqrt((y - yc)*(y - yc) + (x - xc)*(x - xc)) -   r;

 Xm(i) = x0 + Dr(i-1)*cos(theta0) + r * (cos(theta) - cos(theta0));
 Ym(i) = y0 + Dr(i-1)*sin(theta0) + r * (sin(theta) - sin(theta0));
 Zm(i) = z0 + r * tanLamda * (deltaTheta);
 Dz(i) = z - Zm(i) ;
 Zm(i) = Zm(i) + Dz(i);
 %% Prediction



 Ak(:,i) = [Dr(i) theta Dz(i) r  tanLamda]';
 F = [cos(theta - theta0) R0(i)*sin(theta - theta0)  cos(theta - theta0) 0 0;
  	  (-1/R0(i))*sin(theta - theta0) (R0(i-1)/R0(i))*cos(theta-theta0) (-1/R0(i))*sin(theta - theta0) 0 0;
         0 0 1 0 0;
		 (r/R0(i))*tanLamda*sin(theta-theta0) r*tanLamda*(1-(R0(i-1)/R0(i))*cos(theta-theta0)) (r/R0(i))*tanLamda*sin(theta - theta0)  1 -r*deltaTheta;
		 0 0 0 0 1]';

 Fk(1:5,(i-1)*5+1:i*5) = F;
 Ck(1:5,(i-1)*5+1:i*5) = F*Ck1*F';

 %% Filtering 

 Hk = [cos(theta0) sin(theta0) 0;
	  (Dr(i-1)*sin(theta0)-r*sin(theta0)) (r*cos(theta0) - Dr(i-1)*cos(theta0)) 0;
      0 0 1;
      cos(theta)-cos(theta0) sin(theta)-sin(theta0) 0;
	  0 0 r*deltaTheta]';

 Kk = (Ck(1:5,(i-1)*5+1:i*5)^-1 + Hk'*Gk*Hk)^-1 * Hk' * Gk;
 Ck(1:5,(i-1)*5+1:i*5) = (eye(5) - Kk*Hk) * Ck(1:5,(i-1)*5+1:i*5);
 
 Ak(:,i) = Ak(:,i) + Kk * [(X(i) - Xm(i)) (Y(i) - Ym(i)) (Z(i) - Zm(i))]';
 r = Ak(4,i);
 tanLamda = Ak(5,i);

 x0 = x;
 y0 = y;
 z0 = z;
 theta0  = theta;   
endfor

newAk = Ak;
newXm = Xm;
newYm = Ym;
newZm = Zm;

% smoothing
for i=length(X)-1:-1:1
  % update covariance
  NCk = Ck(1:5,(i-1)*5+1:i*5) * Fk(1:5,(i-1)*5+1:i*5)'* Ck(1:5,i*5+1:(i+1)*5);
  
  % update Ak
  newAk(:,i) = Ak(:,i) + NCk * (newAk(:,i+1) - Ak(:,i+1));
endfor

% give back estimate
x0 =M(1,1);
y0 =M(1,2);
z0 =M(1,3);
theta0 = atan2(y0 -yc  ,x0 - xc);
if (theta0 > pi)
	theta0 = theta0 - 2*pi;
endif
if (theta0 < -pi)
	theta0 = 2*pi + theta0;
endif

for i=1:length(X)
 x = M(i,1);
 y = M(i,2);
 z = M(i,3);

 theta = atan2(M(i,2) - yc,M(i,1)-xc);
 % normalize Theta  
 if (theta  > pi)
	theta  = theta  - 2*pi;
 endif
 if (Theta(i) < -pi)
	theta  = theta  + Theta(i);
 endif


 deltaTheta = theta - theta0;
 % normalize deltaTheta  
 if (deltaTheta > pi)
	deltaTheta = deltaTheta- 2*pi;
 endif
 if (deltaTheta < -pi)
	deltaTheta = 2*pi + deltaTheta;
 endif

 Xm(i) = xc + (newAk(4,i) + newAk(1,i)) * cos(newAk(2,i));
 Ym(i) = yc + (newAk(4,i) + newAk(1,i)) * sin(newAk(2,i));
 Zm(i) = z;
 theta0 = theta;

 x0 = x;
 y0 = y;
 z0 = z;

 residue(i) = sqrt((Xm(i) - X(i))^2 + (Ym(i) - Y(i))^2 + (Zm(i) - Z(i))^2);       
endfor

P  = [X' Y' Z'];
Pm = [Xm' Ym' Zm'];
endfunction

