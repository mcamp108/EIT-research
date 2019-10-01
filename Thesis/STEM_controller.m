import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;
% C = get(0, 'PointerLocation');
% screenSize = get(0, 'screensize');
if ~exist('launched', 'var')
    launch_stem();
    launched= 1;
end % end if
n= 2; correct= 10; len= 20;
[nBackSet, answers]= build_nBack_set(n, correct, len);
set_up_screen();

for i= 1:len
    start_capture();
    stimulus= text(0.4, 0.5, nBackSet(i), 'fontsize', 200,'color', 'b');
    pause(1);
    delete(stimulus);
    pause(1);
    stop_capture(answer(i));
end % end for


function launch_stem()
    cd 'C:\Users\Mark\Desktop\STEM';
    system('run.bat');
    cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\Data\Headcap';
    mkdir(date);
end % end function

function start_capture()
    % move mouse to record button (when stem fullscreen)
    mouse.mouseMove(1720, 1080-170);
    left_mouse_click();
end % end function

function stop_capture(isCorrect)
    % called with isCorrect= 1 if stimulus should have elicitied user
    % response
    
    % move mouse to record button (when stem fullscreen)
    mouse.mouseMove(1720, 1080-170);
    left_mouse_click();
    % save file
    if isCorrect
        % name file with should have been positive response
    else
        % name file normally
    end
    
end % end function

function left_mouse_click()
    mouse.mousePress(InputEvent.BUTTON1_MASK);
    mouse.mouseRelease(InputEvent.BUTTON1_MASK);
end % end function

function [nBackSet, answers]= build_nBack_set(n, correct, len)
    % creates an n-back working memory test set from user parameters
    % n: the sequence length of correct answers
    % correct: the number of correct responses
    % length: the total length of the set
    letters= ['A' 'B' 'C' 'D' 'E' 'F' 'G' 'H' 'I' 'J' 'K' 'L' 'M'... 
              'N' 'O' 'P' 'Q' 'R' 'S' 'T' 'U' 'V' 'W' 'X' 'Y' 'Z'];
    % 'length' quantity of random numbers between 1 and 26
    randLetterIdx= randi([1 26], 1, len);
    nBackSet= letters(randLetterIdx);
    % pick a random place for 'correct' number of correct responses to be
    CorrectIdxStart= sort(randi([1 len-n], 1, correct));
    % change the letters 'n' positions away to be the same as those at
    % randCorrectIdx
    CorrectIdxEnd= CorrectIdxStart+ n;
    nBackSet(CorrectIdxEnd)= nBackSet(CorrectIdxStart);
    answers= zeros(1, len);
    answers(CorrectIdxEnd)= 1;
    % check that these are the only correct answers
    dontChange= [CorrectIdxStart, CorrectIdxEnd];
    checknBackSet= 1;
    while checknBackSet== 1
        checknBackSet= 0; % will terminate loop if no corrections are made
        for i= 1:len-n
            if sum(dontChange==i)==0
                if nBackSet(i)== nBackSet(i+n)
                    nBackSet(i)= letters(randi([1 26]));
                    checknBackSet= 1;
                end % end if
            end % end if
        end % end for
    end % end while
    
end % end function

function set_up_screen()
    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'Color', 'white');
    set(gca,'Color','white');
    set(gca,'XColor','white');
    set(gca,'YColor','white');
end % end function
