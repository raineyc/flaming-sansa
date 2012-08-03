function data = scale_to_uint8(data)

% find the mininum height value in the data set;
min_z = min(min(data));

% remove the DC offset in the data
data= data - min_z;

% find the new maximum height;
max_z = max(max(data));

% find the divisor increment
 increment = max_z/256;

%quantize the data 

data = uint8( data/ increment );

return

