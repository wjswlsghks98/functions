function pos_n = pos_update(pos, vel, acc, dt)
%POS_UPDATE : update position in the N frame 
%             (local cartesian frame,assumed inertial)

pos_n = pos + vel*dt + 0.5*acc*dt^2;

end

