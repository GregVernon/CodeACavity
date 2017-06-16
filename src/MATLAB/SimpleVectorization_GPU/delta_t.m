function dt = delta_t(VELOCITY, cREYNOLDS, dt_max, h)

vMax = max(max(VELOCITY));
crMax = max(max(cREYNOLDS));


dt1 = (h^2) / (4*vMax);
dt2 = (2*vMax) / crMax;

dt = min([dt1 dt2 dt_max]);