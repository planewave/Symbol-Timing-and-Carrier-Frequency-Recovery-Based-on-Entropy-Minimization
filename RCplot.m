% to derive eye diagram entropy
close
clear
t= -1:0.01:1;
y1 = sinc(t).*cos(pi*t)./(1-(2*t).^2);
y2 = sinc(t);
y3 = cos(pi*t);
% y3 = y1.*y2;
% plot(t, y1, t, y2, t, y3);
a = 0.4;
y4 = 1 + (4 - (2*pi^2)/3)*t.^2 + 2/15*(120 - 20*pi^2 + pi^4).*t.^4;
% y4 = 1 + (-1/2*(pi^2 - 8)*a^2 - 1/6)*t.^2 + 1/120*(5*(384 - 48*pi^2 + pi^4)*a^4 + 10*(pi^2 - 8)*a^2 + 1)*t.^4 ;
plot(t, y4)