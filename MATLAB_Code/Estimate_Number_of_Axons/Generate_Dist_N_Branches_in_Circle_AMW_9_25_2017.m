% This script places lines (which represent climbing fibers modeled as
% lines of brances) at random angles and displacements relative to a 
% circle (which represents the bottom face of the image volume that is
% modeled as a cylinder).

% The number of interbranch spacings that at maximum will fit in the circle
% (across its diameter) is a fixed value that has been calculated using a
% reasonable value for an inter-branch spacing, along with the size of the
% image volume.

% The number of inter-branch spacings that each modeled climbing fiber has
% inside the circle is calculated based on the length of the line relative
% to the diameter of the circle (which is equal to N).

% The number of branches in two extreme cases is then estimated for
% climbing fiber axons given these numbers of inter-branch spacings for
% each axon.

% These numbers are then used to predict, for any number of
% branches (e.g. the number of branches in the P7 image volume), how many
% of those will come from axons forming 1, 2, 3, etc. branches in the image
% volume.

% Alyssa Wilson, April 2017
% Lichtman Lab, Harvard University

clear all
close all

% Create points along the edge of a circle of diameter N
% Compute the diameter of circle in terms of units, where each unit is the
% width of the interbranch spacing
d_spatial = mean([120.0 75.0]); % Approximate diameter of the circular face
% on the cylindrical image volume
s_spatial = 30.0; % Approximate interbranch spacing from (Sugihara, 2005)
N = d_spatial/s_spatial; % Maximum number of inter-branch spacings/circle diameter
R = N/2;

npoints = 1000;
thetas = rand(npoints,1)*2*pi;

% Convert these theta values to x- y- coordinates
xs = R*cos(thetas);
ys = R*sin(thetas);

figure()
plot(xs,ys,'.');

% Pick the first point in the list as a reference
ref_x = xs(1);
ref_y = ys(1);

% Remove that point from the original list of edge points
xs(1) = [];
ys(1) = [];

% Add a bunch of random points inside the circle.
xs_io = (rand(npoints*10,1)-0.5)*N;
ys_io = (rand(npoints*10,1)-0.5)*N;

rs_io = sqrt(xs_io.*xs_io + ys_io.*ys_io);
outside_circle = find(rs_io > R);
xs_io(outside_circle) = [];
ys_io(outside_circle) = [];

% Look at points generated
figure();
plot(xs_io,ys_io,'.');

% The sticks are the line segments from the reference point to the randomly
% generated points in the circle
stick_length = sqrt((ref_x - xs_io).^2 + (ref_y - ys_io).^2);

% Check the histogram of stick lengths (i.e. of number of inter-branch
% spacings per climbing fiber axon) inside the circle
figure()
bin_edges = [0:1:N];
histogram(stick_length);

% Total number of climbing fiber branches found in the P7 image volume
ntot_branches = 64;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two alternative cases where branches fall at different shifts relative
% to the edge of the circle's edge

% First extreme case: the nearest branch for each climbing fiber is 0.5 the
% inter-branch spacing away from the edge of the circle

% Calculate the number of axons with 1, 2, 3, etc. branches based on the
% rule above
n_zeros_alt1  = length(find(stick_length < 0.5));
n_ones_alt1   = length(intersect(find(stick_length < 1.5),find(stick_length >=0.5)));
n_twos_alt1   = length(intersect(find(stick_length < 2.5),find(stick_length >=1.5)));
n_threes_alt1 = length(intersect(find(stick_length < 3.5),find(stick_length >=2.5)));
n_fours_alt1  = length(intersect(find(stick_length < 4.5),find(stick_length >=3.5)));
n_fives_alt1  = length(find(stick_length >=4.5));

% Calculate the number of branches coming from axons that form 1, 2, 3,
% etc. branches in the volume
n_br_1_alt_1  = n_ones_alt1 * 1;
n_br_2_alt_1  = n_twos_alt1 * 2;
n_br_3_alt_1  = n_threes_alt1 * 3;
n_br_4_alt_1  = n_fours_alt1 * 4;
n_br_5_alt_1  = n_fives_alt1 * 5;

% Calculate the total number of branches in the modeled image volume
n_br_tot_alt_1 = n_br_1_alt_1 + n_br_2_alt_1 + n_br_3_alt_1 + n_br_4_alt_1 + n_br_5_alt_1;

% Calculate the fractions of all branches that come from axons with 1, 2,
% 3, etc. branches in the image volume.
frac_br_1_alt_1 = n_br_1_alt_1/n_br_tot_alt_1
frac_br_2_alt_1 = n_br_2_alt_1/n_br_tot_alt_1
frac_br_3_alt_1 = n_br_3_alt_1/n_br_tot_alt_1
frac_br_4_alt_1 = n_br_4_alt_1/n_br_tot_alt_1
frac_br_5_alt_1 = n_br_5_alt_1/n_br_tot_alt_1

% Calculate the expected number of branches in the P7 volume that come from
% axons with 1, 2, 3, etc. branches in the image volume
n_act_br_1_alt_1 = frac_br_1_alt_1 * ntot_branches;
n_act_br_2_alt_1 = frac_br_2_alt_1 * ntot_branches;
n_act_br_3_alt_1 = frac_br_3_alt_1 * ntot_branches;
n_act_br_4_alt_1 = frac_br_4_alt_1 * ntot_branches;
n_act_br_5_alt_1 = frac_br_5_alt_1 * ntot_branches;

% Convert the expected numbers of branches above into expected numbers of
% axons with 1, 2, 3, etc. branches in the P7 image volume
n_act_ax_1_alt_1 = n_act_br_1_alt_1/1
n_act_ax_2_alt_1 = n_act_br_2_alt_1/2
n_act_ax_3_alt_1 = n_act_br_3_alt_1/3
n_act_ax_4_alt_1 = n_act_br_4_alt_1/4
n_act_ax_5_alt_1 = n_act_br_5_alt_1/5


%%%%%%%%
% Second extreme case: for every climbing fiber, at least one branch is
% located at the edge of the circle

% Calculate the number of axons with 1, 2, 3, etc. branches based on the
% rule above
n_ones_alt_2   = length(find(stick_length < 1));
n_twos_alt_2   = length(intersect(find(stick_length < 2),find(stick_length >=1)));
n_threes_alt_2 = length(intersect(find(stick_length < 3),find(stick_length >=2)));
n_fours_alt_2  = length(intersect(find(stick_length < 4),find(stick_length >=3)));
n_fives_alt_2  = length(find(stick_length >=4));

% Calculate the number of branches coming from axons that form 1, 2, 3,
% etc. branches in the volume
n_br_1_alt_2  = n_ones_alt_2 * 1;
n_br_2_alt_2  = n_twos_alt_2 * 2;
n_br_3_alt_2  = n_threes_alt_2 * 3;
n_br_4_alt_2  = n_fours_alt_2 * 4;
n_br_5_alt_2  = n_fives_alt_2 * 5;

% Calculate the total number of branches in the modeled image volume
n_br_tot_alt_2 = n_br_1_alt_2 + n_br_2_alt_2 + n_br_3_alt_2 + n_br_4_alt_2 + n_br_5_alt_2;

% Calculate the fractions of all branches that come from axons with 1, 2,
% 3, etc. branches in the image volume.
frac_br_1_alt_2 = n_br_1_alt_2/n_br_tot_alt_2
frac_br_2_alt_2 = n_br_2_alt_2/n_br_tot_alt_2
frac_br_3_alt_2 = n_br_3_alt_2/n_br_tot_alt_2
frac_br_4_alt_2 = n_br_4_alt_2/n_br_tot_alt_2
frac_br_5_alt_2 = n_br_5_alt_2/n_br_tot_alt_2

% Calculate the expected number of branches in the P7 volume that come from
% axons with 1, 2, 3, etc. branches in the image volume
n_act_br_1_alt_2 = frac_br_1_alt_2 * ntot_branches;
n_act_br_2_alt_2 = frac_br_2_alt_2 * ntot_branches;
n_act_br_3_alt_2 = frac_br_3_alt_2 * ntot_branches;
n_act_br_4_alt_2 = frac_br_4_alt_2 * ntot_branches;
n_act_br_5_alt_2 = frac_br_5_alt_2 * ntot_branches;

% Convert the expected numbers of branches above into expected numbers of
% axons with 1, 2, 3, etc. branches in the P7 image volume
n_act_ax_1_alt_2 = n_act_br_1_alt_2/1
n_act_ax_2_alt_2 = n_act_br_2_alt_2/2
n_act_ax_3_alt_2 = n_act_br_3_alt_2/3
n_act_ax_4_alt_2 = n_act_br_4_alt_2/4
n_act_ax_5_alt_2 = n_act_br_5_alt_2/5

