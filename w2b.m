function [xb, yb] = w2b(Xt, Xo, heading)
    % World to Body Frame Transformation
    % Xt = [xt, yt] (Target Point coordinates in the Global Frame)
    % Xo = [xv, yv] (Vehicle coordinates in the Global Frame)
    xt = Xt(1); yt = Xt(2);
    xv = Xo(1); yv = Xo(2);
    angle = heading - atan2(yt - yv, xt - xv);
    l = sqrt((xt - xv).^2 + (yt - yv).^2);
    xb = l*cos(angle); yb = l*sin(angle);
end