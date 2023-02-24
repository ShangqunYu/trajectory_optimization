function base = addTerrainToModel(terrain)

%% Grid for parameterized terrain
x_samples = 250;
y_samples = 250;

x_range = [-2.0 5.0];
y_range = [-2.0 2.0];

x = linspace(x_range(1), x_range(2), x_samples);
y = linspace(y_range(1), y_range(2), y_samples);

groundV = zeros(x_samples*y_samples,3);
for iy = 1:y_samples
    for ix = 1:x_samples
        groundV(y_samples*(iy-1)+ix,:) = [x(ix) y(iy) getGroundInfo(terrain, [x(ix) y(iy)])];
    end
end

N_triangles = 2*(x_samples-1)*(y_samples-1);
groundT = zeros(N_triangles, 3);
tri_cnt = 1;
for iy = 1:(y_samples-1)
    for ix = 1:(x_samples-1)
        %groundT(tri_cnt,:) = [y_samples*(iy-1) + ix, y_samples*iy + ix, y_samples*(iy-1) + ix + 1]';
        %groundT(tri_cnt+1,:) = [y_samples*(iy-1) + ix + 1, y_samples*iy + ix, y_samples*iy + ix + 1]';
        groundT(tri_cnt,:) = [y_samples*(iy-1) + ix, y_samples*iy + ix + 1, y_samples*iy + ix]';
        groundT(tri_cnt+1,:) = [y_samples*(iy-1) + ix , y_samples*(iy-1) + ix + 1, y_samples*iy + ix + 1]';
        tri_cnt = tri_cnt + 2;
    end
end

%% Determine Range for tiles
b_valley = false;
x_valley_start = x_range(2);
x_valley_end = x_range(1);
valley_threshold = -0.025;

for ix = 1:x_samples
    
    if ~b_valley && (getGroundInfo(terrain, [x(ix) 0]) < valley_threshold)
        b_valley = true;
        x_valley_start = x(ix);
    end
    
    if b_valley && (getGroundInfo(terrain, [x(ix) 0]) > valley_threshold)
        x_valley_end = x(ix-1);
        break
    end
    
end

%% Create Base
base = {'colour',[0.3 0.3 0.3],...
    'vertices', groundV, ...
    'triangles',groundT,...
    'tiles',[x_range(1) x_valley_start;y_range(1) y_range(2);repmat(min(groundV(:,3))+0.001,1,2)], 0.4,...
    'tiles',[x_valley_end x_range(2);y_range(1) y_range(2);repmat(min(groundV(:,3))+0.001,1,2)], 0.4};
