function [xw, yw] = b2w(Xt, Xo, heading)
    % Body to World Frame Transformation
    % Xt = [xt, yt] (Target Point coordinates in the Body Frame)
    % Xo = [xv, yv] (Vehicle coordinates in the Global Frame)
    xt = Xt(1); yt = Xt(2);
    xv = Xo(1); yv = Xo(2);
    xw = xv + xt * cos(heading) + yt * sin(heading);
    yw = yv + xt * sin(heading) - yt * cos(heading);
end