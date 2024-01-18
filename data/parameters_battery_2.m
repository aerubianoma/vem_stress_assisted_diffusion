function [lambda,mu,m0,m1,phi0,phi1,theta] = parameters_battery_2()
    
    E = 10.e6;
    v = 0.3;
    lambda = (E*v)/((1+v)*(1-2*v));
    mu = 2*(E/(2*(1+v)));

    %% Ricardo's paper

    %m0 = 5.e-1;
    %m1 = 0;
    %phi0 = 1;
    %phi1 = 1;
    %theta = 1.e-3;

    %% Taralova's thesis 

    m0 = 1; %1.e-14;%
    m1 = 0; % 0; % m1 = 0; % For constant M %
    phi0 = ((2*mu+3*lambda)/3)*3.497e-6;
    phi1 = 0; % don't needed in this example %
    theta = 1;