function [rates, err_upp, err_down] = Linear_mass_rates(j,TP1,t)
    ind1 = find(TP1.datenumb<t(j),1,'last');
    ind2 = find(TP1.datenumb>t(j),1,'first');
    day_diff = TP1.datenumb(ind2) - TP1.datenumb(ind1);
    dm = TP1.mass(ind1)-TP1.mass(ind2);
    dt = TP1.datenumb(ind1)-TP1.datenumb(ind2);
    rates = dm/dt;
    %just do the average of two dates (initial and latest)
    iu1 = TP1.error_upp(ind1)/day_diff; id1 = TP1.error_down(ind1)/day_diff;
    iu2 = TP1.error_upp(ind2)/day_diff; id2 = TP1.error_down(ind2)/day_diff;
    err_upp = sqrt(((iu1^2)/4) + ((iu2^2)/4));
    err_down = sqrt(((id1^2)/4) + ((id2^2)/4));
end