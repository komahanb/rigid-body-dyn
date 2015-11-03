function [CBI] = rot(thetain)

    theta  = thetain(3); % rotation in z-plane alone is considered
    
    %C3
         C(1,1) = cos(theta);
         C(1,2) = sin(theta);
         C(1,3) = 0.0d0;
         
         C(2,1) = -sin(theta);
         C(2,2) = cos(theta);
         C(2,3) = 0.0d0;
         
         C(3,1) = 0.0d0;
         C(3,2) = 0.0d0;
         C(3,3) = 1.0d0;
    
    CBI = C;
      
%     C(1,1) = 1.0d0;
%     C(1,2) = 0.0d0;
%     C(1,3) = 0.0d0;
%     
%     C(2,1) = 0.0d0;
%     C(2,2) = cos(theta);
%     C(2,3) = sin(theta);
%     
%     C(3,1) = 0.0d0;
%     C(3,2) = -sin(theta);
%     C(3,3) = cos(theta);
    
%     %C2
%     C(2,1,1) = cos(theta);
%     C(2,1,2) = 0.0d0;
%     C(2,1,3) = -sin(theta);
%     
%     C(2,2,1) = 0.0d0;
%     C(2,2,2) = 1.0d0;
%     C(2,2,3) = 0.0d0;
%     
%     C(2,3,1) = sin(theta);
%     C(2,3,2) = 0.0d0;
%     C(2,3,3) = cos(theta);
%     
%     %C3
%     C(3,1,1) = cos(theta);
%     C(3,1,2) = sin(theta);
%     C(3,1,3) = 0.0d0;
%     
%     C(3,2,1) = -sin(theta);
%     C(3,2,2) = cos(theta);
%     C(3,2,3) = 0.0d0;
%     
%     C(3,3,1) = 0.0d0;
%     C(3,3,2) = 0.0d0;
%     C(3,3,3) = 1.0d0;
%     
             
% 
%     CBI(1,1) =  cos(theta(2))*cos(theta(3)) + sin(theta(1))*sin(theta(2))*sin(theta(3));
%     CBI(2,1) = -cos(theta(2))*sin(theta(3)) + sin(theta(1))*sin(theta(2))*cos(theta(3));
%     CBI(3,1) =  cos(theta(1))*sin(theta(2));
% 
%     CBI(1,2) =  cos(theta(1))*sin(theta(3));
%     CBI(2,2) =  cos(theta(1))*cos(theta(3));
%     CBI(3,2) = -sin(theta(1));
% 
%     CBI(1,3) = -sin(theta(2))*cos(theta(3)) + sin(theta(1))*cos(theta(2))*sin(theta(3));
%     CBI(2,3) =  sin(theta(2))*sin(theta(3)) + sin(theta(1))*cos(theta(2))*cos(theta(3));
%     CBI(3,3) =  cos(theta(1))*cos(theta(2));
%     
%     %%  Eq. 19 Hughes
%     CBI(1,1) =  cos(theta(2))*cos(theta(3));
%     CBI(1,2) =  cos(theta(2))*sin(theta(3));
%     CBI(1,3) =  -sin(theta(2));
% 
%     CBI(2,1) =  sin(theta(1))*sin(theta(2))*cos(theta(3)) - cos(theta(1))*sin(theta(3));
%     CBI(2,2) =  sin(theta(1))*sin(theta(2))*sin(theta(3)) + cos(theta(1))*cos(theta(3));
%     CBI(2,3) =  sin(theta(1))*cos(theta(2));
% 
%     CBI(3,1) =  cos(theta(1))*sin(theta(2))*cos(theta(3)) + sin(theta(1))*sin(theta(3));
%     CBI(3,2) =  cos(theta(1))*sin(theta(2))*sin(theta(3)) - sin(theta(1))*cos(theta(3));
%     CBI(3,3) =  cos(theta(1))*cos(theta(2));
    
   %page 20
   % write( new rotation matrix and smatrix and perhaps sdot) 
   % check the code execution
   % check the approximation now
   % try to extract the s_approx
   % use S instead of s_approx
end