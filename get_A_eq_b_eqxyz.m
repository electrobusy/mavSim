function [A_x,b_x]=get_A_eq_b_eqxyz(keyframes,t,x,dx,ddx,dddx,ddddx,variablename)
m=size(keyframes,2);
A_p_test=zeros(4*(m-1),length(x)*(m-1));
A_p_eq=[];
b_p_eq=[];
n=length(x)-1;

if strcmp(variablename,'x')
    j=1;
elseif strcmp(variablename,'y')
    j=2;
elseif strcmp(variablename,'z')
    j=3;
elseif strcmp(variablename,'psi')
    j=4;
end
for i=1:(m-1) %loop over number of sections between keyframes
    
   %% Equality constraints (x=keyframe, [dx=0, ddx=0 (for begin and end caps)] 
   A_p_i=[
       polyval_terms(x,t(i)); 
       polyval_terms(x,t(i+1));];
   b_p_eq=[b_p_eq;
       keyframes(j,i);
       keyframes(j,i+1)];     
   if i==1 %if the first or last polynomial add dx and ddx equality constraints
       A_p_i=[A_p_i;
           [polyval_terms(dx,t(i)),0];
           [polyval_terms(ddx,t(i)),0,0]
           ];
       b_p_eq=[b_p_eq;0;0];
   end
   if i==(m-1)
       A_p_i=[A_p_i;
           [polyval_terms(dx,t(m)),0];
           [polyval_terms(ddx,t(m)),0,0]
           ];
       b_p_eq=[b_p_eq;0;0];
   end
    A_p_eq=blkdiag(A_p_eq,A_p_i); 
end

%% Continuity constraints at the inner waypoints
if strcmp(variablename,'psi')
    A_p_cont=zeros(3*(m-2),2*(n+1)+(m-3)*(n+1));
    b_p_cont=zeros(size(A_p_cont,1),1);
    for i =1:(m-2) % loop over inner keyframes for continuity constraints 
        A_p_cont((1+(i-1)*3):(3+(i-1)*3),(1+(i-1)*(n+1)):(2*(n+1)+(i-1)*(n+1)))=[ %replace 5 with the number of continuity constraints if it dynamically changes
        polyval_terms(x,t(i+1)), -polyval_terms(x,t(i+1));
        polyval_terms(dx,t(i+1)) 0, -polyval_terms(dx,t(i+1)) 0;
        polyval_terms(ddx,t(i+1)) 0 0, -polyval_terms(ddx,t(i+1)) 0 0];
    end
    
else
    A_p_cont=zeros(5*(m-2),2*(n+1)+(m-3)*(n+1));
    b_p_cont=zeros(size(A_p_cont,1),1);
    for i =1:(m-2) % loop over inner keyframes for continuity constraints 
        A_p_cont((1+(i-1)*5):(5+(i-1)*5),(1+(i-1)*(n+1)):(2*(n+1)+(i-1)*(n+1)))=[ %replace 5 with the number of continuity constraints if it dynamically changes
        polyval_terms(x,t(i+1)), -polyval_terms(x,t(i+1));
        polyval_terms(dx,t(i+1)) 0, -polyval_terms(dx,t(i+1)) 0;
        polyval_terms(ddx,t(i+1)) 0 0, -polyval_terms(ddx,t(i+1)) 0 0;
        polyval_terms(dddx,t(i+1)) 0 0 0, -polyval_terms(dddx,t(i+1)) 0 0 0;
        polyval_terms(ddddx,t(i+1)) 0 0 0 0, -polyval_terms(ddddx,t(i+1)) 0 0 0 0];
    end
end


A_x=[A_p_eq;A_p_cont];
b_x=[b_p_eq;b_p_cont];
end