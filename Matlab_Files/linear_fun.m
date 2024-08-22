%Linearly interp
function rates = linear_fun (j,TP,t)
dm = TP.mass(find(TP.datenumb<t(j),1,'last'))-TP.mass(find(TP.datenumb>t(j),1,'first'));
dt = TP.datenumb(find(TP.datenumb<t(j),1,'last'))-TP.datenumb(find(TP.datenumb>t(j),1,'first'));
rates = dm/dt;
end

