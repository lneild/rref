% Homework Program 5
%
% Name: Neild, Lainey
% Date: October 22, 2021

function [M, piv] = reduce(M)

% initialize pivots to no pivot cols so far - each val = col number that contains a piv
piv = []; % impossible to preallocate space - unclear number of pivs

% not returned by function but used to get to rref from echelon form 
piv_row_col = []; 

% initialize the starting point to find piv 1
cur_row = 1;
cur_col = 1;

% this loop = zeros under all pivs

while (cur_row <= size(M,1) && cur_col <= size(M,2))
    
    %redefine our matrix of M to the submatrix
    sub_mat = M(cur_row:end, cur_col:end);
    
    % find the max elem in each col of sub matrix
    [max_per_col, pos_of_max] = max(abs(sub_mat), [], 1);
    
    % find row associated with first nonzero elem in max_per_col
    % tells you piv col of row - if no piv in col will move to next col
    col_sub_mat = find(max_per_col,1); % col with first pivot
    if isempty(col_sub_mat)
       break;
    end
    piv_row_sub_mat = pos_of_max(col_sub_mat); % row of submatrix with max
     
    % find col of M from col of sub_mat
    M_col = col_sub_mat + cur_col - 1;
    
    % redfine submatrix - in case row of all 0s
    sub_mat = M(cur_row:end, M_col:end);
    
    piv = [piv, M_col]; % add this col number to the piv vector
    
    % update piv_row_col
    piv_row_col = [piv_row_col; cur_row M_col]; % in reverse order
    
    % swaping rows in submatrix and in matrix M
    
    % swap row in sub matrix
    sub_mat = exchange(sub_mat, piv_row_sub_mat);
    
    % swap row in M matrix - M_row is the row you are swappign with cur_row
    M_row = cur_row + piv_row_sub_mat - 1;
    M = exchange(M, cur_row, M_row);
    
    row_of_piv = cur_row; % update row_of_piv now that you have switched it
    pivot_val = sub_mat(1, 1); % piv val - used to find val to multiply
   
    % values under piv - will 0 out
    entries_under_piv = sub_mat(2:end, 1);
    
    % number of nonzero entries under piv 
    entries_below = nnz(entries_under_piv);
    
    % row of entries currently zeroing out
    find(entries_under_piv,1);
    temp_row = 1 + find(entries_under_piv,1);
    
    while ((entries_below > 0)) % while still have vals left to 0 out
        
        % what multiply to get to zero out - add by -1 = subtract
        factor = sub_mat(temp_row,1)/pivot_val*-1;
        
        %zero out in sub matrix
        sub_mat = add(sub_mat, factor, 1, temp_row);
        
        % zero out in M
        M = add(M, factor, row_of_piv, (temp_row + cur_row - 1));
        
        % we have to assign it to the value of zero -  round off error
        % if value is -0.0000 will set entries remaining to 1 - error!
        sub_mat(temp_row,1) = 0; % change to 0 in submatrix
        M((temp_row + cur_row - 1), M_col) = 0; % change to 0 in M
        
        
        % update entries under piv to exclude the one just zeroed
        entries_under_piv = sub_mat(temp_row:end, 1);
        
        %update number of rows left to zero out
        entries_below = nnz(entries_under_piv);
        
        % update your temp_row to be next nonzero entry under piv
        temp_row = temp_row + 1;
        
        % if next row = 0, skip over it
        while (entries_below > 0 && (sub_mat(temp_row,1) == 0))
            temp_row=temp_row+1;
        end
    end
    
    % move to next row/col as starting point to find next piv
    cur_row = cur_row + 1;
    cur_col = M_col + 1; % in case skipped a col
   
    
end

% this loop zeros above all pivs - goes through each of col with a piv
% start at 2 - nothing to zero above the first piv
for k = 2:(1+size(piv_row_col,1))
   
    % current row and col of piv val - need to zero above M(r,c) in col c
    % subtract 1 from k because we started k at 2
    r = piv_row_col((k-1),1);
    c = piv_row_col((k-1),2);
    p_val = M(r,c);
    
    % make pivot 1 by multiplying the entire row by 1/p_val
    M(r,:) = M(r,:)*(1/p_val);
    
    % this zeros out each row in each col goes from 1 above piv row to 1
    for current_r = (r-1):-1:1
        
        if (M(current_r,c) == 0)
            continue; % if this val is already zeroed out, no operation
        end
        
        % what multiply to get to zero out - add by -1 = subtract
        f = M(current_r,c)*-1;
        
        % zero out in M
        M = add(M, f, r, current_r);
        
        % we have to assign it to the value of zero -  round off error
        % if value is -0.0000 will set entries remaining to 1 - error!
        M(current_r, c) = 0; % change to 0 in M
      
    end
    
end

end % end of func reduce


function M = exchange (M, row1, row2)
% exchanges rows in matrix M

% if both empty set both to default
if (~exist('row1', 'var') || isempty(row1)) && (~exist('row2', 'var') || isempty(row2))
    row1 = 1; 
    row2 = 2; 
    
% if enetered a var for r1 and not for r2, exchange row 1 and r1
elseif exist('row1', 'var') && ~exist('row2', 'var')
    row2 = row1;
    row1 = 1;

end

if row1 ~= round(row1) % check if  r1 is int
    error('This is not an integer');

elseif row1 < 1 || row1 > size(M, 1) % if int, check if valid
    error('This is not a valid input');
end

if row2 ~= round(row2) % check if r2 is int
    error('This is not an integer');

elseif row2 < 1 || row2 > size(M, 1) % if int, check if valid
    error('This is not a valid input');
    
end

% tell user what rows are being exchanged
fprintf('This function will exchange row %i and row %i.\n', row1, row2);

% exchange the rows
M([row1 row2], :) = M([row2 row1], :);

end % end of func exchange


function M = mult (M, d, row)
% multiplies row r by nonzero value d in matrix M, returns updated M

% check if d = 0
if (d == 0)
    error('You cannot enter zero for this value');
end

% check if didn't enter val set to default
if ~exist('row', 'var') || isempty(row)
    row = 1; % set to default val
end

if row ~= round(row) % is row an int
    error('You have to enter an integer for the row number');
    
elseif row < 1 || row > size(M, 1) % if did enter val, is it valid
    error('This is not valid input');
end


% tell user which row is getting multiplied and by what
fprintf('Row %i is getting multiplied by %.4f\n', row, d);

% multiply row by d in M 
M(row, :) = M(row, :) *d;

end % end of func mult

function M = add(M, r, row1, row2)
% adds r times one row and adds it to another row in a matrix
% replaces row 2 with r times row 1

% if enetered r1 and not for r2, multiply r by row entered, add to row 1 
if exist('row1', 'var') && ~exist('row2', 'var')
    row2 = 1;
end

if row1 ~= round(row1) % check if  r1 is int
    error('This is not an integer');
elseif row1 < 1 || row1 > size(M, 1) % check if valid row number
    error('This is not valid input');
end

if row2 ~= round(row2) % check if  r2 is int
    error('This is not an integer');
elseif row2 < 1 || row1 > size(M, 1) % check if valid row number
    error('This is not valid input');
end

% tell user what is happening
fprintf('Row %i is being multiplied by %.4f and added to row %i\n', ...
    row1, r, row2);

M(row2, :) = M(row2, :) + r*M(row1, :);

end % end of func add


% outputs:
% 
% 1. A = randi([-5 5], 4, 10); [R,piv] = reduce(A)
% This function will exchange row 1 and row 3.
% This function will exchange row 1 and row 3.
% Row 1 is being multiplied by 0.3333 and added to row 2
% Row 1 is being multiplied by 0.3333 and added to row 2
% Row 1 is being multiplied by 0.3333 and added to row 3
% Row 1 is being multiplied by 0.3333 and added to row 3
% Row 1 is being multiplied by -1.0000 and added to row 4
% Row 1 is being multiplied by -1.0000 and added to row 4
% This function will exchange row 1 and row 2.
% This function will exchange row 2 and row 3.
% Row 1 is being multiplied by -0.1000 and added to row 2
% Row 2 is being multiplied by -0.1000 and added to row 3
% Row 1 is being multiplied by 0.9000 and added to row 3
% Row 2 is being multiplied by 0.9000 and added to row 4
% This function will exchange row 1 and row 2.
% This function will exchange row 3 and row 4.
% Row 1 is being multiplied by -0.4231 and added to row 2
% Row 3 is being multiplied by -0.4231 and added to row 4
% This function will exchange row 1 and row 1.
% This function will exchange row 4 and row 4.
% Row 2 is being multiplied by 0.3333 and added to row 1
% Row 3 is being multiplied by 0.4000 and added to row 2
% Row 3 is being multiplied by 0.8000 and added to row 1
% Row 4 is being multiplied by -0.8846 and added to row 3
% Row 4 is being multiplied by -0.1538 and added to row 2
% Row 4 is being multiplied by 0.6923 and added to row 1
% 
% R =
% 
%     1.0000         0         0         0    0.1623    0.9215   -0.7749   -0.0890    0.8272    0.1675
%          0    1.0000         0         0   -1.5916   -0.6492   -1.2723    0.1309   -1.6283    0.5183
%          0         0    1.0000         0   -0.1518    0.2670    0.4346   -1.4974   -0.6126   -0.7696
%          0         0         0    1.0000    0.3455    0.7199   -0.7304   -0.3508    0.0838    0.1309
% 
% 
% piv =
% 
%      1     2     3     4
%
% 2. A = randi([-5 5], 10, 4); [R, piv] = reduce(A)
% This function will exchange row 1 and row 7.
% This function will exchange row 1 and row 7.
% Row 1 is being multiplied by 0.2000 and added to row 2
% Row 1 is being multiplied by 0.2000 and added to row 2
% Row 1 is being multiplied by 0.6000 and added to row 3
% Row 1 is being multiplied by 0.6000 and added to row 3
% Row 1 is being multiplied by 0.2000 and added to row 4
% Row 1 is being multiplied by 0.2000 and added to row 4
% Row 1 is being multiplied by 0.8000 and added to row 5
% Row 1 is being multiplied by 0.8000 and added to row 5
% Row 1 is being multiplied by 0.8000 and added to row 6
% Row 1 is being multiplied by 0.8000 and added to row 6
% Row 1 is being multiplied by -0.6000 and added to row 7
% Row 1 is being multiplied by -0.6000 and added to row 7
% Row 1 is being multiplied by -1.0000 and added to row 8
% Row 1 is being multiplied by -1.0000 and added to row 8
% Row 1 is being multiplied by -0.2000 and added to row 9
% Row 1 is being multiplied by -0.2000 and added to row 9
% Row 1 is being multiplied by 1.0000 and added to row 10
% Row 1 is being multiplied by 1.0000 and added to row 10
% This function will exchange row 1 and row 2.
% This function will exchange row 2 and row 3.
% Row 1 is being multiplied by 0.3077 and added to row 2
% Row 2 is being multiplied by 0.3077 and added to row 3
% Row 1 is being multiplied by 0.8846 and added to row 3
% Row 2 is being multiplied by 0.8846 and added to row 4
% Row 1 is being multiplied by 0.6538 and added to row 4
% Row 2 is being multiplied by 0.6538 and added to row 5
% Row 1 is being multiplied by 0.4615 and added to row 5
% Row 2 is being multiplied by 0.4615 and added to row 6
% Row 1 is being multiplied by 0.8077 and added to row 6
% Row 2 is being multiplied by 0.8077 and added to row 7
% Row 1 is being multiplied by -0.1923 and added to row 7
% Row 2 is being multiplied by -0.1923 and added to row 8
% Row 1 is being multiplied by -0.3077 and added to row 8
% Row 2 is being multiplied by -0.3077 and added to row 9
% Row 1 is being multiplied by -0.1923 and added to row 9
% Row 2 is being multiplied by -0.1923 and added to row 10
% This function will exchange row 1 and row 8.
% This function will exchange row 3 and row 10.
% Row 1 is being multiplied by -0.1667 and added to row 2
% Row 3 is being multiplied by -0.1667 and added to row 4
% Row 1 is being multiplied by 0.4286 and added to row 3
% Row 3 is being multiplied by 0.4286 and added to row 5
% Row 1 is being multiplied by -0.4167 and added to row 4
% Row 3 is being multiplied by -0.4167 and added to row 6
% Row 1 is being multiplied by 0.5476 and added to row 5
% Row 3 is being multiplied by 0.5476 and added to row 7
% Row 1 is being multiplied by 0.2381 and added to row 6
% Row 3 is being multiplied by 0.2381 and added to row 8
% Row 1 is being multiplied by 0.3810 and added to row 7
% Row 3 is being multiplied by 0.3810 and added to row 9
% Row 1 is being multiplied by -0.2262 and added to row 8
% Row 3 is being multiplied by -0.2262 and added to row 10
% This function will exchange row 1 and row 4.
% This function will exchange row 4 and row 7.
% Row 1 is being multiplied by -0.0430 and added to row 2
% Row 4 is being multiplied by -0.0430 and added to row 5
% Row 1 is being multiplied by 0.4892 and added to row 3
% Row 4 is being multiplied by 0.4892 and added to row 6
% Row 1 is being multiplied by 0.2258 and added to row 4
% Row 4 is being multiplied by 0.2258 and added to row 7
% Row 1 is being multiplied by -0.1075 and added to row 5
% Row 4 is being multiplied by -0.1075 and added to row 8
% Row 1 is being multiplied by -0.6237 and added to row 6
% Row 4 is being multiplied by -0.6237 and added to row 9
% Row 1 is being multiplied by -0.3495 and added to row 7
% Row 4 is being multiplied by -0.3495 and added to row 10
% Row 2 is being multiplied by -0.4000 and added to row 1
% Row 3 is being multiplied by -0.4615 and added to row 2
% Row 3 is being multiplied by 0.3846 and added to row 1
% Row 4 is being multiplied by 0.4643 and added to row 3
% Row 4 is being multiplied by -0.2143 and added to row 2
% Row 4 is being multiplied by 0.1786 and added to row 1
% 
% R =
% 
%      1     0     0     0
%      0     1     0     0
%      0     0     1     0
%      0     0     0     1
%      0     0     0     0
%      0     0     0     0
%      0     0     0     0
%      0     0     0     0
%      0     0     0     0
%      0     0     0     0
% 
% 
% piv =
% 
%      1     2     3     4
% 
% 3. A = [1 2 0 0 0;0 0 0 2 3;0 0 2 1 3]; [R, piv] = reduce(A)
% This function will exchange row 1 and row 1.
% This function will exchange row 1 and row 1.
% This function will exchange row 1 and row 2.
% This function will exchange row 2 and row 3.
% This function will exchange row 1 and row 1.
% This function will exchange row 3 and row 3.
% Row 3 is being multiplied by -0.5000 and added to row 2
% 
% R =
% 
%     1.0000    2.0000         0         0         0
%          0         0    1.0000         0    0.7500
%          0         0         0    1.0000    1.5000
% 
% 
% piv =
% 
%      1     3     4
%
% 4. a = [1 2 4 3 5]; A = [a;a;a;a]; [R, piv] = reduce(A)
% 
% This function will exchange row 1 and row 1.
% This function will exchange row 1 and row 1.
% Row 1 is being multiplied by -1.0000 and added to row 2
% Row 1 is being multiplied by -1.0000 and added to row 2
% Row 1 is being multiplied by -1.0000 and added to row 3
% Row 1 is being multiplied by -1.0000 and added to row 3
% Row 1 is being multiplied by -1.0000 and added to row 4
% Row 1 is being multiplied by -1.0000 and added to row 4
% 
% R =
% 
%      1     2     4     3     5
%      0     0     0     0     0
%      0     0     0     0     0
%      0     0     0     0     0
% 
% 
% piv =
% 
%      1
