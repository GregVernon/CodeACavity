function dt = delta_t(VELOCITY, cREYNOLDS, dt_max, h, NY, NX)

vMax = 0;
crMax = 0;
for jj = 1:NY
    for ii = 1:NX
        if VELOCITY(jj,ii) > vMax
            vMax = VELOCITY(jj,ii);
        end
        if cREYNOLDS(jj,ii) > crMax
            crMax = cREYNOLDS(jj,ii);
        end
    end
end

dt1 = (h^2) / (4*vMax);
dt2 = (2*vMax) / crMax;

dt = dt1;
if dt2 < dt
    dt = dt2;
elseif dt_max < dt
    dt = dt_max;
end
    

