function [ out ] = isProcessOrMovieList( x )
%isProcessOrMovieData True if input is a Process or MovieData instance

out = isProcessOrMovieObject(x, 'MovieList');

end
