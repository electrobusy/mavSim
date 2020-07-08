function [t,poly] =  discretize_poly(coeffs,T,dt)
    nr_sections=size(coeffs,3);
    poly_order=size(coeffs,1);
    nr_states=size(coeffs,2);
    nr_derivative=size(coeffs,4)-1;

    t=nan((max(diff(T)))/dt+1,poly_order,nr_sections); %each section may have a different length, the ends are padded with nans
    % build time array to multiply coefficients with
    poly=nan(size(t,1),nr_states,nr_sections,nr_derivative+1);
    for s = 1:nr_sections
        ts=T(s):dt:T(s+1);
        for n=1:size(t,2)
            t(1:length(ts),n,s)=ts.^(n-1);
        end
        for d=1:nr_derivative+1
            poly(:,:,s,d)=t(:,:,s)*coeffs(:,:,s,d);
        end
    end
    t=t(:,2,:);
    end