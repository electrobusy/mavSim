function [t,poly] =  discretize_poly(coeffs,T,data)
    dt=data.dt;
    nr_sections=size(coeffs,3);
    poly_order=size(coeffs,1);
    nr_states=size(coeffs,2);
    nr_derivative=size(coeffs,4)-1;
    tl=[];
    dT=diff(T);
    t=nan((max(diff(T)))/dt,poly_order,nr_sections); %each section may have a different length, the ends are padded with nans
    % build time array to multiply coefficients with
%     poly=nan(size(t,1),nr_states,nr_sections,nr_derivative+1);
    poly=nan((T(end)-T(1))/dt,nr_states,nr_derivative+1);
    k=1;
    for s = 1:nr_sections
        if isfield(data,'anoosh')
            ts=0:dt:dT(s)-dt;
            ts2=T(s):dt:T(s+1)-dt;
        else
            ts=T(s):dt:T(s+1)-dt;
        end
        for n=1:size(t,2)
            t(1:length(ts),n,s)=ts.^(n-1);
        end
        for d=1:nr_derivative+1
            vals=t(1:length(ts),:,s)*coeffs(:,:,s,d);
            poly(k:size(vals,1)+k-1,:,d)=vals;
        end
        k=k+size(ts,2);
        
        if isfield(data,'anoosh')
           tl=[tl,ts2];
        else
            tl=[tl,ts];
        end
    end
 
    
    t=tl;%reshape(t(:,2,:),size(poly,1),1);
    end