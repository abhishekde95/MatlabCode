% J Patrick Weller
% 2/4/12
% PBio 545

% Created   1/4/12  JPW

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem Set #1

% Question set #1:

% (1) Why does the identity transform a vector back into itself?  Pick
%       some examples and convince yourself this is true. 

% The identity matrix transforms a vector back into itself because it
% reproduces each element of the vector independently - that is, when taking
% the dot product, the identity matrix multiplies some element of the vector
% by one, and adds to it every other element of the vector multiplied by
% zero, thereby preserving its initial identity.  Since each element of the
% vector undergoes this same transformation, the final vector is the same as
% the original vector.

% [1 0 0    [x     [x+0+0  
%  0 1 0  *  y  =   0+y+0 
%  0 0 1]    z]     0+0+z]

IDmat = eye(3,3)
randmat = unifrnd(1,10,3,1)

transformedMat = IDmat*randmat


% (2) Confirm that the u we get above is right by the brute force approach
%		to the original coupled linear equations.

% 2x + 3y - 3z = 1
% -x + y = 2
% 3x - y + z = -1

% Multiply the third equation by 3 and add it to the first equation:
% 2x + 3y - 3z + 9x - 3y + 3z = -3 + 1
% 11x = -2
% x = -2/11 (~ -0.1818)

% Plug this value for x into the second equation:
% -(-2/11) + y = 2
% y = 20/11 (~ 1.8182)

% Plug x and y into the third equation:
% 3(-2/11) - (20/11) + z = -1
% -6/11 - 20/11 + 1 = -z
% z = 15/11 (~ 1.3636)

% Thus, we confirm by the brute force approach that the vector u, which we
% obtained using the inverse M matrix, is correct.


% (3) Why does checking that M*u = v tell us we found the right solution
%       u?

% Crudely speaking, checking that M*u = v tells us that we have found the
% correct solution because we were looking for a vector (u) that is
% transformed into vector v when the operator M is applied to it.  Since
% our vector u does become vector v when M is applied to it, we know that
% we have found the correct solution.

% Further, we know that we have found the correct solution because M is an
% invertible matrix - that is, one that yields a unique 1-1 mapping of one
% vector into another, and does not, alternatively, transform several
% different vectors (u', u'', u''', etc) into the same vector v.



%% Question set #2:


% (1) Why do rotation matrices take the form above (e.g. the 2-d rotation matrix)?  
%	  Draw a picture and make use of the definitions of sine and cosine from trigonometry.

%1. Rotation matrices take the form [cos? -sin?; sin? cos?] because this is
%the general solution for rotating any vector about any angle ? while
%preserving the length of that vector.  Here is a proof:

theta = unifrnd(0,360)
gamma = unifrnd(0,360)
x = unifrnd(1,10)
y = unifrnd(1,10)


% To be proven:  norm([cos? -sin?; sin? cos?] * [x; y] ) = norm([x;y])

(cos(theta)*x - sin(theta)*y).^2 + (sin(theta)*x + cos(theta)*y).^2
cos(theta).^2*x.^2 - 2*x*y*cos(theta)*sin(theta) + sin(theta).^2*y.^2 + sin(theta).^2*x.^2 + 2*x*y*cos(theta)*sin(theta) + cos(theta).^2*y.^2
(cos(theta).^2*x.^2 + sin(theta).^2*x.^2) + (sin(theta).^2*y.^2 + cos(theta).^2*y.^2)
x.^2*(cos(theta).^2 + sin(theta).^2) + y.^2*(sin(theta).^2 + cos(theta).^2)
	% since cos(?).^2 + sin(?).^2 = 1
%norm([cos(theta) -sin(theta); sin(theta) cos(theta)] * [x; y] ) = x.^2+ y.^2
norm([cos(theta) -sin(theta); sin(theta) cos(theta)] * [x; y])
sqrt(x.^2 + y.^2)

% Thus, norm([cos? -sin?; sin? cos?] * [x; y] ) = norm([x;y]), which was to be proven.

% QED



% (2) What would the matrix look like for a rotation about the x-axis of 30 degrees
%     followed by a rotation about the z-axis of 45 degrees?

% Rotation about the x-axis of 30 degrees:

xRotMat = [1	 0	 	0
        0 	cos(30) -sin(30)
        0 	sin(30)  cos(30)]

% Rotation about the z-axis of 45 degrees:

zRotMat = [cos(45) -sin(45)	0
            sin(45) cos(45)	0
            0       0		1]

%Since these are both linear transformations, both matrices can be combined into one.  Order does matter here.

rotscalemat =    [cos(theta) -sin(theta)*cos(gamma) sin(theta)*sin(gamma)
                 sin(theta) cos(theta)*cos(gamma) -cos(theta)*sin(gamma)
                    0            sin(gamma)             cos(gamma)]


% (3) When a rotation transformation and a scaling transformation are both applied to a
%	  vector, does the order make a difference?  Can you show this?  

% When a rotation transformation and a scaling transformation are both
% applied to a vector,  the order does make a difference.

theta = unifrnd(0,360);
scalemat = [.5 0; 0 2];
rotmat = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% When the scaling matrix is applied to the rotation matrix, we get:
scale2Rot = scalemat * rotmat

% When, however, the rotation matrix is applied to the scaling matrix, we get:
rot2Scale = rotmat * scalemat

% Thus, order does matter when applying linear transformations.


% (4) Why do we use the inverse of M above as a scaling matrix?  Why does the inverse
%	  take the form it does?

% We apply the inverse of M and apply it to the matrix 'Dists' because M
% was the matrix by which we transformed the orignal draws from a normal
% Gaussian distribution.  Thus, in order to transform Dists back into its
% original form (that is, in which the distance of a point from the origin
% is equal to its standard deviation), the inverse of M needs to be applied
% to it.


%% Question Set #3

% (1) Give an example of a matrix of rank 1.  What are the column space and null spaces for this 
%	    matrix?

rank1Mat = [1 0 0
            0 0 0
            0 0 0];
rank(rank1Mat)
        
colSpace = orth(rank1Mat)
        
nullSpace = null(rank1Mat)



% (2) How many linearly independent vectors are represented in a matrix of rank 5?

% In a matrix of rank 5, 5 linearly independent vectors are represented.



