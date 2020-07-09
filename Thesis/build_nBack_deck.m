function [n_back_deck, answers]= build_nBack_deck(n, correct, len)

% -------------------------------------------------------------------------
% DESCRIPTION:
%   [n_back_deck, answers]= build_nBack_deck(n, correct, len)
%   
%   creates an n-back working memory test set from user parameters
% -------------------------------------------------------------------------
% PARAMETERS:
%   n: 
%       the sequence length of correct answers
%   correct: 
%       the number of correct responses
%   length: 
%       the total length of the deck
% -------------------------------------------------------------------------   
% RETURNS:
%   n_back_deck:
%
%   answers:
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------

letters= ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M'... 
          'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];

% 'length' quantity of random numbers between 1 and 26
randLetterIdx   = randi([1 26], 1, len);
n_back_deck     = letters(randLetterIdx);

% pick a random place for 'correct' number of correct responses to be
CorrectIdxStart = sort(randi([1 len-n], 1, correct));

% change the letters 'n' positions away to be the same as those at
% randCorrectIdx
CorrectIdxEnd   = CorrectIdxStart + n;
n_back_deck(CorrectIdxEnd)= n_back_deck(CorrectIdxStart);
answers         = zeros(1, len);
answers(CorrectIdxEnd)= 1;

% check that these are the only correct answers
dontChange      = [CorrectIdxStart, CorrectIdxEnd];
checknBackSet   = 1;
while checknBackSet == 1
    checknBackSet = 0; % will terminate loop if no corrections are made
    for i = 1:len-n
        if sum(dontChange==i) ==0
            if n_back_deck(i) == n_back_deck(i+n)
                n_back_deck(i) = letters(randi([1 26]));
                checknBackSet = 1;
            end % end if
        end % end if
    end % end for
end % end while

end % end function