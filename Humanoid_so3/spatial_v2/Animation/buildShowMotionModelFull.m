function model = buildShowMotionModelFull(params)

switch params.model
    case 'quad3D'
        model = buildShowMotionModelMC3DFull(params);
    case 'humanoid3D'
        model = buildShowMotionModelH3DFull(params);
    case 'humanoidNoArms3D'
        model = buildShowMotionModelHNA3DFull(params);
end

% Ground
if nargin == 2
    model.appearance.base = {'tiles', [-10 20;-6 6;ground_height ground_height], 0.4};
else
    model.appearance.base = {'tiles', [-10 20;-6 6;0 0], 0.4};
end

% Camera
model.camera.body = 6;
model.camera.direction = [0.05 -0.6 0.2];
model.camera.zoom = 0.65;