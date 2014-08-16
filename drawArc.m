function drawArc(c, radius, start_angle, stop_angle)

angle = linspace(start_angle, stop_angle, 1000);

X = c.x + radius.*cos(angle);
Y = c.y + radius.*sin(angle);

plot(X, Y);


p1.x = X(1);
p1.y = Y(1);

p2.x = X(length(X));
p2.y = Y(length(Y));

drawLine(c, p1);
drawLine(c, p2);


end