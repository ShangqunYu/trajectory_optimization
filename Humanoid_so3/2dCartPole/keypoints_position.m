function out1 = keypoints_position(in1,in2)
%KEYPOINTS_POSITION
%    OUT1 = KEYPOINTS_POSITION(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    24-Nov-2022 18:09:48

lbar = in2(3,:);
th = in1(2,:);
x = in1(1,:);
out1 = reshape([x+lbar.*sin(th),-lbar.*cos(th),x,0.0],[2,2]);