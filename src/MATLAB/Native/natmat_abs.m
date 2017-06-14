function v = natmat_abs(v)
% https://stackoverflow.com/questions/9772348/
v = v * ((v>0) - (v<0));
% if v < 0
%     v = -v;
% end