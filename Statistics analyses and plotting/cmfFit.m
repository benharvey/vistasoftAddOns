function B=cmfFit(ii,x,y,ve)
x = x(ii); y = y(ii); ve = ve(ii);
B = fminsearch(@(z) mycmffit(z,x,y,ve),[0.05;0.2]);
return

function e=mycmffit(z,x,y,ve)
e=sum(ve.*(y-(1./(z(1).*x+z(2)))).^2)./sum(ve);
return
