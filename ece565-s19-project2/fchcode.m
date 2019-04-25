function c = fchcode(b, CONN) 
%FCHCODE Computes the Freeman chain code of a boundary. 

%C = FCHCODE(B, CONN) produces the same outputs as above, but 
% with the code connectivity specified in CONN. CONN can be 8 for 
% an 8-connected chain code, or CONN can be 4 for a 4-connected 
% chain code. Specifying CONN=4 is valid only if the input 
% sequence, B, contains transitions with values 0, 2, 4, and 6, 
% exclusively. % % C = FHCODE(B, CONN, DIR) produces the same outputs as above, but, 
% in addition, the desired code direction is specified. Values for % DIR can be: 
% 
% 'same' Same as the order of the sequence of points in b. 
% This is the default. 
% 
% 'reverse' Outputs the code in the direction opposite to the 
% direction of the points in B. The starting point % for each DIR is the same. 
% 
% The elements of B are assumed to correspond to a 1-pixel-thick, 
% fully-connected, closed boundary. B cannot contain duplicate % coordinate pairs, except in the first and last positions, which 
% is a common feature of boundary tracing programs. 
% 
% FREEMAN CHAIN CODE REPRESENTATION 
% The table on the left shows the 8-connected Freeman chain codes 
% corresponding to allowed deltax, deltay pairs. An 8-chain is 
% converted to a 4-chain if (1) if conn = 4; and (2) only 
% transitions 0, 2, 4, and 6 occur in the 8-code. Note that % dividing 0, 2, 4, and 6 by 2 produce the 4-code. 
% % ----------------------- ---------------- 
% deltax | deltay | 8-code corresp 4-code 
% ----------------------- ---------------- 
% 0 1 0 0 
% -1 1 1 
% -1 0 2 1 
% -1 -1 3 
% 0 -1 4 2 
% 1 -1 5 
% 1 0 6 3 
% 1 1 7 
% ----------------------- ---------------- 
% 
% Preliminaries. 
if nargin == 1
    dir = 'same'; 
    conn = 8; 
elseif nargin == 2 
    dir = 'same'; 
elseif nargin == 3 
    % Nothing to do here. 
else
    error('Incorrect number of inputs.')
end
[np, nc] = size(b); 
if np < nc 
    error('B must be of size np-by-2.'); 
end
% Some boundary tracing programs, such as boundaries.m, output a
% sequence in which the coordinates of the first and last points are 
% the same. If this is the case, eliminate the last point. 
if isequal(b(1, :), b(np, :)) 
    np = np - 1; 
    b = b(1:np, :); 
end
% Build the code table using the single indices from the formula 
% for z given above: 
C(11)=0;
C(7)=1;
C(6)=2;
C(5)=3;
C(9)=4;
C(13)=5;
C(14)=6;
C(15)=7; 
% End of Preliminaries. 
% Begin processing. 
x0 = b(1, 1); 
y0 = b(1, 2); 
c.x0y0 = [x0, y0]; 
% Make sure the coordinates are organized sequentially: 
% Get the deltax and deltay between successive points in b. The 
% last row of a is the first row of b. 
a = circshift(b, [-1, 0]); 
% DEL = a - b is an nr-by-2 matrix in which the rows contain the 
% deltax and deltay between successive points in b. The two 
% components in the kth row of matrix DEL are deltax and deltay 
% between point (xk, yk) and (xk+1, yk+1). The last row of DEL 
% contains the deltax and deltay between (xnr, ynr) and (x1, y1), 
% (i.e., between the last and first points in b). 
DEL = a - b; 
% If the abs value of either (or both) components of a pair 
% (deltax, deltay) is greater than 1, then by definition the curve 
% is broken (or the points are out of order), and the program 
% terminates. 
if any(abs(DEL(:, 1)) > 1) | any(abs(DEL(:, 2)) > 1); 
    error('The input curve is broken or points are out of order.') 
end 
% Create a single index vector using the formula described above. 
z = 4*(DEL(:, 1) + 2) + (DEL(:, 2) + 2); 
% Use the index to map into the table. The following are 
% the Freeman 8-chain codes, organized in a 1-by-np array. 
fcc = C(z) ;
% Check if direction of code sequence needs to be reversed. 
if strcmp(dir, 'reverse')
fcc = coderev(fcc); 
% See below for function coderev. 
end 
% If 4-connectivity is specified, check that all components 
% of fcc are 0, 2, 4, or 6. 
if CONN == 4 
    val = find(fcc == 1 | fcc == 3 | fcc == 5 | fcc ==7 );
    if isempty(val) 
        fcc = fcc./2; 
    else
        warning('The specified 4-connected code cannot be satisfied.')
    end
end
% Freeman chain code for structure output. 
c.fcc = fcc;

% Obtain the first difference of fcc. 
c.diff = firstdiff(fcc,CONN);
% See below for function codediff. 
% Obtain code of the integer of minimum magnitude. 
c.mm = min_mag(fcc);
% See below for function minmag. 
% Obtain the first difference of fcc 
c.diffmm = firstdiff(c.mm, CONN);
disp(c.diff)
end