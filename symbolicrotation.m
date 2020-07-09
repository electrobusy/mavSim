function R=symbolicrotation(angles,axorder,orientation,demo)
%% this function produces the rotation matrix from a sequence of unit rotations
% The orientation states whether it concerns a lefthanded or righthanded axis system
% demo is for symbolic representation to check with paper
% angles should be like [ax, ay, az] regardless of what axorder is
% pad angles with zeros if using less than 3 angles 
% currently only works with 3 rotations or less 
% angles
% demo=1
% orientation='righthanded';
% axorder=[3 2 1]; %note that this order represents the order of angle rotations is the reverse of the matrix multiplication order
if demo
syms cphi sphi ctheta stheta cpsi spsi
sphi=sign(angles(1))*sphi;
stheta=sign(angles(2))*stheta;
spsi=sign(angles(3))*spsi;
else
    cphi=cos(angles(1));
    sphi=sin(angles(1));
    ctheta=cos(angles(2));
    stheta=sin(angles(2));
    cpsi=cos(angles(3));
    spsi=sin(angles(3));
end
if strcmp(orientation,'righthanded')
   if ismember(1,axorder)
      Rx=[1, 0 ,0;
          0, cphi, -sphi;
          0, sphi, cphi];  
   else 
       Rx=eye(3);      
   end
   if ismember(2,axorder)
       Ry=[ctheta, 0, stheta; 
           0,  1 , 0;
      
           -stheta, 0 , ctheta];
   else
       Ry=eye(3);
   end
   if ismember(3,axorder)
       Rz=[cpsi, -spsi, 0;
           spsi, cpsi, 0;
           0 ,  0 , 1]; 
   else
       Rz=eye(3);
   end
   
elseif strcmp(orientation,'lefthanded') %left handed inverts the sign of the sinus terms
    if ismember(1,axorder)
      Rx=[1, 0 ,0;
          0, cphi, sphi;
          0, -sphi, cphi];  
   else 
       Rx=eye(3);      
   end
   if ismember(2,axorder)
       Ry=[ctheta, 0, -stheta; 
           0,  1 , 0;      
           stheta, 0 , ctheta];
   else
       Ry=eye(3);
   end
   if ismember(3,axorder)
       Rz=[cpsi, spsi, 0;
           -spsi, cpsi, 0;
           0 ,  0 , 1]; 
   else
       Rz=eye(3);
   end
else
    error("Orientation not recognized \n Choose 'lefthanded' or 'righthanded'");
end

unitrotmats=cat(3,Rx,Ry,Rz);
R=eye(3);
for i=flip(axorder) %should order be flipped? (according to flight dynamics course it should, but some papers don't seems to do it)
    R=R*unitrotmats(:,:,i);
end
% R
% end