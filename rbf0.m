function phi = rbf0(c,t,ti)

r = t - ti;
R = abs(r);

h = 1/c*R;

F = sqrt(h.^2+1);

%Multiquadric
phi = F + R.^5;

end