function [CBI] = rot(theta)

    CBI(1,1) =  cos(theta(2))*cos(theta(3)) + sin(theta(1))*sin(theta(2))*sin(theta(3));
    CBI(2,1) = -cos(theta(2))*sin(theta(3)) + sin(theta(1))*sin(theta(2))*cos(theta(3));
    CBI(3,1) =  cos(theta(1))*sin(theta(2));

    CBI(1,2) =  cos(theta(1))*sin(theta(3));
    CBI(2,2) =  cos(theta(1))*cos(theta(3));
    CBI(3,2) = -sin(theta(1));

    CBI(1,3) = -sin(theta(2))*cos(theta(3)) + sin(theta(1))*cos(theta(2))*sin(theta(3));
    CBI(2,3) =  sin(theta(2))*sin(theta(3)) + sin(theta(1))*cos(theta(2))*cos(theta(3));
    CBI(3,3) =  cos(theta(1))*cos(theta(2));
   %page 20
   % write( new rotation matrix and smatrix and perhaps sdot) 
   % check the code execution
   % check the approximation now
   % try to extract the s_approx
   % use S instead of s_approx
end