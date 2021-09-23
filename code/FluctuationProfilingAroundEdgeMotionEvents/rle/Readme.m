% Unzip all files in Matlab current directory .
%If you have your vector for example:
%X=[56 56 56 2 28 28 255 255 9 9 9];%This is the which we want to cmpress it.
%Now we can call to our function :
%rle(X)
%ans=
%    56  3  2  1  28  2  255  2  9  3
%Here X is the Input and ans is the Output, in output 56 is the value of
%pixel and 3 is the successive répitition.
%We can use this function only with values of vetors, do not use it with
%characters.
%rle.m : Run Length Encoding function
%irle.m : Inverse Run Length Encoding function
%K=[85 3 41 1 96 4 210 1];
%irle(K)
%ans=
%    85  85  85  41  96  96  96  96  210 