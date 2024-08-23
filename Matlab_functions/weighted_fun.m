%Weighted Interp
function rates = weighted_fun (terms,TP,date_sum,~)
        clear rates avg_mass_rate sum_avg_mass_rate;
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
        last_date_diff = last_day-TP.datenumb(terms(l));
%         else
%             mass_rate(l+1) = 0;
%             last_date_diff = 0;
%        end
        first_date_diff = TP.datenumb(terms(1))-first_day;
%       Differencing from the first day of the month
        if first_date_diff == 0
            avg_mass_rate = mass_rate(1)*1;
        else
            avg_mass_rate = mass_rate(1)*first_date_diff;
        end

        for i = 1:length(terms)-1
            avg_mass_rate = avg_mass_rate + date_change(i)*mass_rate(i+1);
        end
        if last_date_diff == 0
%           last_date_change = TP.datenumb(terms(l)+1) - TP.datenumb(terms(l));
            sum_avg_mass_rate = avg_mass_rate + (mass_rate(l+1)* 1);
        else
            sum_avg_mass_rate = avg_mass_rate + (mass_rate(l+1)* last_date_diff);
        end
%         date_sum = first_date_diff + last_date_diff;
%         for j = 1:length(date_change)
%             date_sum = date_sum + date_change(j);
%         end
        rates = sum_avg_mass_rate/date_sum;
        
end


