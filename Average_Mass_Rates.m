%Average mass rates
function [rates,max_rate,min_rate] = Average_Mass_Rates(terms,TP,~)
        clear rates avg_mass_rate sum_avg_mass_rate mass_rate;
        l = length(terms);
        first_day = datenum(dateshift(TP.dates(terms(1)),'start','month'));
        last_day = datenum (dateshift(TP.dates(terms(1)),'end','month'));
        first_mass_change = TP.mass(terms(1)) - TP.mass(terms(1)-1);
        first_date_change = TP.datenumb(terms(1)) - TP.datenumb(terms(1)-1);
        mass_change = zeros(l,1);
        date_change = zeros(l,1);
        mass_rate = zeros(l+1,1);
        mass_rate(1,1) = first_mass_change/first_date_change;
        if l>1
            for k = 1:l-1
                mass_change(k) = TP.mass(terms(k+1)) - TP.mass(terms(k));
                date_change(k) = TP.datenumb(terms(k+1)) - TP.datenumb(terms(k));
                mass_rate(k+1) = mass_change(k)/date_change(k);
            end
        end
        last_mass_change = TP.mass(terms(l)+1) - TP.mass(terms(l));
        last_date_change = TP.datenumb(terms(l)+1) - TP.datenumb(terms(l));
        mass_rate(l+1) = last_mass_change/last_date_change;
        rates = mean(mass_rate);
        max_rate = max(mass_rate);
        min_rate = min(mass_rate);
end



