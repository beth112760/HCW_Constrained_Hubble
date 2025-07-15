function phid = rbf1(c,t,ti)


r = t - ti;
R = abs(r);

h = 1/c*R;

%Multiquadric
phid = (h./c./sqrt(h.^2+1) + 5*R.^4).*sign(r);



end