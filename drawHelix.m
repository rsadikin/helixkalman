
clear all;
close all;

%A =csvread("1541HelixZPos.csv");
A =csvread("1540Helix.csv"); 
%A =csvread("2039Helix.csv"); %Bad Track
%A =csvread("2590Helix.csv"); %Bad Track
%A =csvread("2798Helix.csv"); %Bad Track
%A =csvread("687Helix.csv");  %Bad Track
%A =csvread("761Helix.csv");  %Bad Track

X= A(:,1)';
Y= A(:,2)';
Z= A(:,3)';

figure
subplot(1,2,1) 
plot3(Z,X,Y,'r*')
title('3D View')

subplot(1,2,2)
plot(X,Y,'b.')
title('2D Front View')

