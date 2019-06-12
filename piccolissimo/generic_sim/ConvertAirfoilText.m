function [ out_mat ] = ConvertAirfoilText( in_string )
%ConvertAirfoilText( in_string ) takes in a string of values from the
%sandia Aerodynamic characteristics of seven symmetrical airfoil sections
%paper and tries to convert it to real data.

in_string = strrep(in_string, '—','-');
in_string = strrep(in_string, ',','');

parsed_nums = textscan(in_string,'%f');

[out_mat,padded] = vec2mat(parsed_nums{1},4);
if padded
    warning('Matrix may not be aligned. Check input string.');
end

out_mat(:,1) = [];
end

