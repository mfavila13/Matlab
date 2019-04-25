function c = fchcode(b, CONN) 
%dfine indexing
C(11)=0;
C(7)=1;
C(6)=2;
C(5)=3;
C(9)=4;
C(13)=5;
C(14)=6;
C(15)=7; 

x0 = b(1, 1); 
y0 = b(1, 2); 
% reference point
c.x0y0 = [x0, y0]; 

a = circshift(b, [-1, 0]); 
DEL = a - b;
z = 4*(DEL(:, 1) + 2) + (DEL(:, 2) + 2);

% use x and y diffference to assign numerical direction
fcc = C(z);

% convert to 4-connected if needed
if CONN == 4 
    val = find(fcc == 1 | fcc == 3 | fcc == 5 | fcc ==7 );
    if isempty(val) 
        fcc = fcc./2; 
    else warning('The specified 4-connected code cannot be satisfied.') 
    end
end
% freeman chain code
c.fcc = fcc;

% first difference of c.fcc
c.diff = firstdiff(fcc,CONN);

% integer of minimum magnitude of c.fcc
c.mm = min_mag(fcc);

% first difference of c.mm 
c.diffmm = firstdiff(c.mm, CONN);bound2im
end