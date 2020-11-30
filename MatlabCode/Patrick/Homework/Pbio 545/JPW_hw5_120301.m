% J Patrick Weller
% PBio 545 Homework Assignment #5
% Created   3/1/12     JPW

clear all
close all

%% Question set #1:

%	(1) Explain the value of the components of the 2 dimensional covariance matrix in 
%		Example 2 above. 

Dist(:, 1) = normrnd(0,0.5,10000,1);
Dist(:, 2) = normrnd(0,1.5,10000,1);

Dist(:, 3) = normrnd(0,0.5,10000,1);			
Dist(:, 4) = normrnd(0,1.5,10000,1) + 1.5 * Dist(:, 3);

figure(1); clf;
subplot(1, 2, 1)
plot(Dist(:, 1), Dist(:, 2), '.');
axis([-6 6 -6 6])
axis square
xlabel('x')
ylabel('y')
title('Uncorrelated')
subplot(1, 2, 2)
plot(Dist(:, 3), Dist(:, 4), '.');
axis([-6 6 -6 6])
axis square
xlabel('x')
ylabel('y')
title('Correlated')

Covar1 = cov(Dist(:,1),Dist(:,2))
Covar2 = cov(Dist(:,3),Dist(:,4))

% The top left and bottom right component in this 2 dimensional covariance
% matrix are the variances of each of the individual distributions (x and
% y).  These values of Covar1 should not surprise us because they are the
% variances that we asked Matlab to create into the distributions (although
% we used standard deviations instead of variance: 0.5 and 1.5).  The upper
% right and lower left values are artifacts of having a finite amount of
% data.  While our two random variables have no actual covariance, any
% finite amount of draws from those distributions will have some near zero
% covariance.
% The variance of x in Covar2 is the same as in Covar1, but the variance of
% y is larger.  This is because we have added an additional amount of
% variance to each of our y components: we have added scaled values from
% our x distribution.  The upper right and lower left components are
% non-zero because we have a real covariance between x and y.
% Additionally, we have a positive value because there is a positive
% covariance between x and y.  That is to say, when x is positive, there is
% a greater-than-chance probability that y is also positive; and when x is
% negative, there is the same probability that y is also negative.  Thus,
% when we multiply each of these x and y values together, then sum, we have
% some positive value.  If there were no covariance, and when x was
% positive, y is just as likely to be positive as negative, these would sum
% to zero (with round-off error).  Note that we also adjust the
% distributions for their absolute positions within the axes - this is so
% we don't mistake two distributions without any covariance that both live
% entirely in the positive domain for having a positive covariance.

%%
%	(2) Which of the off diagonal elements in the covariance matrix for Example 4 above
%		are nonzero?  Why?

clear Dist;
NumDim = 6
for dim = 1:NumDim
	Dist(:, dim) = normrnd(0,dim,50000, 1);			
	if (dim > 1)
		Dist(:, dim) = Dist(:, dim) + 0.5 * Dist(:, dim-1);			
	end
end

cov(Dist)

% All of the off-diagonal elements in this covariance matrix are non-zero.
% This is because each element in a given dimension has some positive
% correlation with each element of the previous (points already drawn from
% the distribution of another dimension).  Thus, neighboring dimensions
% have the highest covariance: x has the greatest covariance with y ((1,2) 
% or (2,1) in the matrix), dimension 6 has the greatest covariance with
% dimension 5 ((6,5) or (5,6) in the matrix), and less covariance with 
% other dimensions the farther away they fall in the matrix.

%%
%	(3) Can you separate inv(Covar^0.5) as above into the product of a rotation and 
%		a scaling matrix?

% Yes, it should be possible.  Both rotation and scaling are linear
% processes, so they are seperable.  

clear Dist

Dist(:, 1) = normrnd(0,0.5,10000,1);			
Dist(:, 2) = normrnd(0,0.5,10000,1) + 1.5 * Dist(:, 1);

Covar = cov(Dist);

ScaleMat = [1/Covar(1,1)^.5 0; 0 1/Covar(2,2)^.5];
ScaleDist = Dist * ScaleMat';
CovSc = cov(ScaleDist);

RotMat = inv(Covar^.5)/ScaleMat;
RotDist =  Dist * RotMat;
CovRot = cov(RotDist);

RotScMat = inv(Covar^.5);
RotScDist = Dist * RotScMat;
CovRS = cov(RotScDist);

figure(3);
subplot(2,2,1)
plot(Dist(:,1),Dist(:,2),'.')
axis([-6 6 -6 6])
axis square
title('Original Distribution')
text(-5,5,['Var X = ' num2str(Covar(1,1))])
text(-5,4,['Var Y = ' num2str(Covar(2,2))])
text(-5,3,['Covariance = ' num2str(Covar(1,2))])
subplot(2,2,2)
plot(RotDist(:,1),RotDist(:,2),'.')
axis([-6 6 -6 6])
axis square
title('Rotated Distribution')
text(-5,5,['Var X = ' num2str(CovRot(1,1))])
text(-5,4,['Var Y = ' num2str(CovRot(2,2))])
text(-5,-5,['Covariance = ' num2str(CovRot(1,2))])
subplot(2,2,3)
plot(ScaleDist(:,1),ScaleDist(:,2),'.')
axis([-6 6 -6 6])
axis square
title('Scaled Distribution')
text(-5,5,['Var X = ' num2str(CovSc(1,1))])
text(-5,4,['Var Y = ' num2str(CovSc(2,2))])
text(-5,3,['Covariance = ' num2str(CovSc(1,2))])
subplot(2,2,4)
plot(RotScDist(:,1),RotScDist(:,2),'.')
axis([-6 6 -6 6])
axis square
title('Rotated and Scaled Distribution')
text(-5,5,['Var X = ' num2str(CovRS(1,1))])
text(-5,4,['Var Y = ' num2str(CovRS(2,2))])
text(-5,-5,['Covariance = ' num2str(CovRS(1,2))])


%% Question set #2:

%	(1) Can you think of another case where dimensional reduction is important?

% Anytime we are trying to find important components in a many-component
% system, dimension reduction will be useful.  The examples in the papers
% we are reading this week are clear examples: we have a cell that has some
% preferred stimulus.  The space of possible preferred stimuli is quite
% large (depending on how many things you vary, such as color, spacial
% extent or location, motion speed, motion direction, etc).  One
% possibility for finding the preferred stimulus (or, more likely, class of
% stimuli) without presenting every possible stimulus is to sample this
% space randomly, then find the directions of greatest variance.  These
% are the dimensions that the cell seems to "care about" the most.

%	(2) How much of the variance of the singles are we capturing with the first
%	    3 components above?  How much more do we get with 10 components?

ThreeDimVar = sum(sum(EigVal(118:120,118:120)));
TotalVar = sum(sum(EigVal));
ThreeDimVar/TotalVar

% With only three dimensions, we account for more than 83% of the total
% variance.

TenDimVar = sum(sum(EigVal(111:120,111:120)))
TenDimVar/TotalVar

% Using ten dimensions, we account for over 93% of the total variance.


%   (3) Why do we subtract the covariance of the failures from that of
%       the singles? 

% We subtract the covariance of the failures from the singles because we
% assume that the latter is still present in the former, even though a
% photon has increased the response of the cell.  Thus, if we did not
% subtract this "dark variance" from the activity of a cell exposed to a
% photon, then we would be mixing the variance of one with the variance of
% another, and not isolating the variance due to the activation by a
% photon.


%   (4) Why don't any of the eigenvectors we recover look like the mean single
%       photon response?

% The eigenvectors don't look like the mean single photon response because
% each of the eigenvectors are measures of variance, not components of any
% of the individual responses.  An eigenvector tells us which direction
% (from the origin) along which a collection of points has the greatest
% variance.  Since that group of points may be indefinitely far from the
% origin, and since thier greatest variance may be in a direction
% orthoganal to the direction of the origin, that eigenvector will bear 
% very little resemblance to the mean of those points.