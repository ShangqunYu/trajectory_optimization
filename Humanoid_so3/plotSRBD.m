function plotSRBD(varargin)
% PLOTCUBE - Display a 3D-cube in the current axes
%
%   PLOTCUBE(EDGES,ORIGIN,ALPHA,COLOR) displays a 3D-cube in the current axes
%   with the following properties:
%   * EDGES : 3-elements vector that defines the length of cube edges
%   * ORIGIN: 3-elements vector that defines the start point of the cube
%   * ALPHA : scalar that defines the transparency of the cube faces (from 0
%             to 1)
%   * COLOR : 3-elements vector that defines the faces color of the cube
%
% Example:
%   >> plotcube([5 5 5],[ 2  2  2],.8,[1 0 0]);
%   >> plotcube([5 5 5],[10 10 10],.8,[0 1 0]);
%   >> plotcube([5 5 5],[20 20 20],.8,[0 0 1]);
% Default input arguments
inArgs = { ...
  [10 56 100] , ... % Default edge sizes (x,y and z)
  [10 10  10] , ... % Default coordinates of the origin point of the cube
  [0 0 0;
  0 0 0;
  0 0 0]     , ... % rotation matrix  
  .4          , ... % Default alpha value for the cube's faces
  [0 0 0.5]       ... % Default Color for the cube
  };
% Replace default input arguments by input values
inArgs(1:nargin) = varargin;
% Create all variables
[edges,origin,rot_m, alpha,clr] = deal(inArgs{:});

rx = [edges(1);0;0]; ry = [0;edges(2);0]; rz = [0;0;edges(3)];

%8 vertices in local frame
ov1 = -1/2 *rx - 1/2 * ry - 1/2 * rz;
ov2 = +1/2 *rx - 1/2 * ry - 1/2 * rz;
ov3 = +1/2 *rx - 1/2 * ry + 1/2 * rz;
ov4 = -1/2 *rx - 1/2 * ry + 1/2 * rz;

ov5 = -1/2 *rx + 1/2 * ry - 1/2 * rz;
ov6 = +1/2 *rx + 1/2 * ry - 1/2 * rz;
ov7 = +1/2 *rx + 1/2 * ry + 1/2 * rz;
ov8 = -1/2 *rx + 1/2 * ry + 1/2 * rz;


%8 vertices in world frame
v1 = origin' + rot_m * ov1;
v2 = origin' + rot_m * ov2;
v3 = origin' + rot_m * ov3;
v4 = origin' + rot_m * ov4;
v5 = origin' + rot_m * ov5;
v6 = origin' + rot_m * ov6;
v7 = origin' + rot_m * ov7;
v8 = origin' + rot_m * ov8;

% six patches
% x, y ,z coordinates for each surface.
patch([v1(1),v2(1),v3(1),v4(1)],[v1(2),v2(2),v3(2),v4(2)],[v1(3),v2(3),v3(3),v4(3)],clr, 'FaceAlpha',alpha);  % right
patch([v1(1),v4(1),v8(1),v5(1)],[v1(2),v4(2),v8(2),v5(2)],[v1(3),v4(3),v8(3),v5(3)],[0.5 0.5 0.5], 'FaceAlpha',0.95);  % back
patch([v5(1),v6(1),v7(1),v8(1)],[v5(2),v6(2),v7(2),v8(2)],[v5(3),v6(3),v7(3),v8(3)],clr, 'FaceAlpha',alpha);  % left
patch([v2(1),v3(1),v7(1),v6(1)],[v2(2),v3(2),v7(2),v6(2)],[v2(3),v3(3),v7(3),v6(3)],clr, 'FaceAlpha',alpha);  % front
patch([v3(1),v4(1),v8(1),v7(1)],[v3(2),v4(2),v8(2),v7(2)],[v3(3),v4(3),v8(3),v7(3)],clr, 'FaceAlpha',alpha);  % top
patch([v1(1),v2(1),v6(1),v5(1)],[v1(2),v2(2),v6(2),v5(2)],[v1(3),v2(3),v6(3),v5(3)],clr, 'FaceAlpha',alpha);  % bottom

view(3);