function [edges_bin, edges_angle, edges_angle_im] = find_canny_edges(data, operator ,threshold)
%% find_canny_edges
% Author: Cameron Rainey
% Contact : raineyc@vt.edu
%
% This function uses canny edge detection as described on 
% http://en.wikipedia.org/wiki/Canny_edge_detector
% 
% USAGE:
% [edges_bin, edges_angle, edges_angle_im] = find_canny_edges(data, method, threshold)
%
% Arguments:
% data - a 2D matrix or image. 
% operator - This specifies the edge detection operator to use.
%            Supported opperators include 'roberts', 'prewitt', 'sobel' and
%            'mat_grad', which uses the matlab gradient function.
% threshold - a threshold value to surpress noise that can appear in the
%             final edged data. (by default the median value of data is used)
%
% Returns:
% edges_bin - a 2D matrix of equal size to the input matrix 'data'.  This
%             matrix contains a binary values, where 0 = no edge, and 1 =
%             edge.
% edges_angle - a 2D matrix of equal size to the input matrix 'data'.  This
%               matrix contains the angle of the gradient (see wikipedia
%               page) rounded to one of the following angles: 0 -> 0
%               degrees, 1 -> 45 degrees, 2-> 90 degrees, 3-> 135 degrees.
%               
% edges_angle_im - returns a 3D matrix, or more precisely, the RGB matrix
%                  of the edge angles where 0 degrees is red, 45 degrees is yellow, 90
%                  degrees is blue, and 135 degrees is green.  This can be
%                  displayed using the 'imshow' command.

if ~exist('operator','var');
    operator = 'mat_grad';
end





operator = lower(operator);

[G, theta] = find_G_theta(data,operator);

if ~exist('threshold','var');
    threshold = median(median(G));
end




G(G<threshold) = 0;


quantized_theta = mod(round(theta/(pi/4)),4);

[edges_bin, edges_angle, edges_angle_im] = find_edges(data, G, quantized_theta);



end

function [G, theta] = find_G_theta(data, operator)

switch operator
    case 'roberts'
        
        % Define Roberts operators
        R1 = [1 0; 0 -1];
        R2 = [0 1; -1 0];
        
        % apply Roberts operators to find Gi and Gj
        Gj = conv2(data,R1,'same');
        Gi = conv2(data,R2,'same');
        
        
    case 'prewitt'
        
        % Define Prewitt operators
        P1 = [-1 0 1; -1 0 1; -1 0 1];
        P2 = [1 1 1; 0 0 0; -1 -1 -1];
        
        
        % apply Prewitt operators to find Gi and Gj
        Gj = conv2(data,P1,'same');
        Gi = conv2(data,P2,'same');
        
        
    case 'sobel'
        % Define Sobel operators
        S1 = [-1 0 1; -2 0 2; -1 0 1];
        S2 = [-1 -2 -1; 0 0 0; 1 2 1];
        
        % apply Sobel operators to find Gi and Gj
        Gj = conv2(data,S1,'same');
        Gi = conv2(data,S2,'same');
        
        
    case 'mat_grad'
        
        % find gradient
        [Gj, Gi] = gradient(data);
        
    otherwise
        warning(strcat('Unrecognised edge detection operator: ', operator, '. Defaulting to mat_grid'));
        
        % find gradient
        [Gj, Gi] = gradient(data);
        
end



% calculate G
G = sqrt(Gj.^2 + Gi.^2);

% Calculate theta values
theta = atan2(Gi,Gj);

% unwrap theta values
% *** This could probably be better accomplished ***
theta(theta < 0) = theta(theta < 0) + 2*pi;

return
end

function [edges_bin, edges_angle, edges_angle_im] = find_edges(data, G, quantized_theta)


 
% Find the size the data matrix
[rows, cols] = size(data);

% allocate the output matricies
edges_bin = false([rows, cols] );             % This matrix is for storing the binary edge values
edges_angle = uint8(zeros([rows, cols] ));    % This matrix stores the angles of the edges
edges_angle_im = uint8(zeros(rows, cols,3 )); % This matrix stores the angles of the edges as a colorcoded image

i_limit = rows;
j_limit = cols;


for i = 1:i_limit
    for j = 1:j_limit
        
        switch quantized_theta(i,j)
            case 0
                
                %if the rounded gradient angle is zero degrees (i.e. the edge is in the
                % north-south direction) the point will be considered to be on the edge if
                % its gradient magnitude is greater than the magnitudes in the west and east
                % directions,
                
                
                j_ind1 = j-1;
                i_ind1 = i;
                j_ind3 = j+1;
                i_ind3 = i;
                
                % Check to see if points are out of bounds.  If points are
                % out of bounds, then set the point values to be NaN
                if j_ind1 < 1, j_ind1 = NaN; end
                if j_ind3 > j_limit, j_ind3 = NaN; end
                
                angle_color = [255,0,0];
                
            case 1
                % if the rounded gradient angle is 135 degrees (i.e. the edge is in the
                % north east-south west direction) the point will be considered to be on
                % the edge if its gradient magnitude is greater than the magnitudes in the
                % north west and south east directions,
                
                j_ind1 = j-1;
                i_ind1 = i-1;
                j_ind3 = j+1;
                i_ind3 = i+1;
                
                
                % Check to see if points are out of bounds.  If points are
                % out of bounds, then set the point values to be NaN
                if j_ind1 < 1, j_ind1 = NaN; end
                if j_ind3 > j_limit, j_ind3 = NaN; end
                if i_ind1 < 1, i_ind1 = NaN; end
                if i_ind3 > i_limit, i_ind3 = NaN; end
                
                angle_color = [255, 255, 0];
                
                
                
                
            case 2
                % if the rounded gradient angle is 90 degrees (i.e. the edge is in the east
                % -west direction) the point will be considered to be on the edge if its
                % gradient magnitude is greater than the magnitudes in the north and south
                % directions,
                
                j_ind1 = j;
                i_ind1 = i-1;
                j_ind3 = j;
                i_ind3 = i+1;
                
                
                % Check to see if points are out of bounds.  If points are
                % out of bounds, then set the point values to be NaN
                if i_ind1 < 1, i_ind1 = NaN; end
                if i_ind3 > i_limit, i_ind3 = NaN; end
                
                angle_color = [0, 0, 255];
                
            case 3
                % if the rounded gradient angle is 45 degrees (i.e. the edge is in the
                % north west-south east direction)the point will be considered to be on the
                % edge if its gradient magnitude is greater than the magnitudes in the north
                % east and south west directions.
                j_ind1 = j-1;
                i_ind1 = i+1;
                j_ind3 = j+1;
                i_ind3 = i-1;
                
                
                % Check to see if points are out of bounds.  If points are
                % out of bounds, then set the point values to be NaN
                if j_ind1 < 1, j_ind1 = NaN; end
                if j_ind3 > j_limit, j_ind3 = NaN; end
                if i_ind1 > i_limit, i_ind1 = NaN; end
                if i_ind3 < 1, i_ind3 = NaN; end
                
                angle_color = [0,255,0];
                
                
                
                
                
                
            otherwise
                warning(['Unexpected theta value: ' num2str(quantized_theta(i,j))]);
         
        end
        
        % Matrix of the points to retrieve
        inds = [i_ind1, j_ind1 ;...
                  i_ind3, j_ind3,];
        
              
        % Determine if any of the points are out of bounds. if so, remove
        % them.
        
        failed_inds = isnan(sum(inds,2));
        
        % Remove the offending points
        inds = inds(~failed_inds,:);
        
        
        
        if ~isempty(inds)
            % determine if this is the maximum point
            % this check requires that it be equal to or greater than one
            % of the two check pixles, however, all three values cannot be
            % the same.
            is_max = ~sum(G(i,j) < G(sub2ind([rows,cols], inds(:,1), inds(:,2)))) && sum(G(i,j) ~= G(sub2ind([rows,cols], inds(:,1), inds(:,2))));
            
            if is_max
                % assign a binary edge
                edges_bin(i,j) = true;
                % assign an edge which marks the angle
                edges_angle(i,j) = quantized_theta(i,j);
                
                edges_angle_im(i,j,:) = angle_color;
                
            end
        end
                
        
    end
end
end







