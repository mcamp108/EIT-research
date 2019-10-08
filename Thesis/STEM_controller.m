import java.awt.Robot;
import java.awt.event.*;
Mouse.mouse = Robot;
Mouse.left_click = InputEvent.BUTTON1_MASK;
Mouse.right_click = InputEvent.BUTTON2_MASK;
% C = get(0, 'PointerLocation');
screenSize = get(0, 'screensize');
if ~exist('launched', 'var')
    launch_stem(Mouse);
    launched= 1;
end % end if
n= 2; correct= 10; len= 20;
[nBackSet, answers]= build_nBack_set(n, correct, len);
set_up_screen();
% save_button= [1722 168];
% confirm_button= [1353 493];
for i= 1:len
    start_capture(Mouse);
    stimulus= text(0.4, 0.5, nBackSet(i), 'fontsize', 200,'color', 'b');
    pause(1);
    delete(stimulus);
    pause(1);
    response= '';
    stop_capture(Mouse, answers(i), response);
    keyboard;
end % end for


function launch_stem(Mouse)
    cd 'C:\Users\Mark\Desktop\STEM';
    system('run.bat');
    cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\Data\Headcap';
    mkdir(date);
    % click sensor quality screen
    Mouse.mouse.mouseMove(232, 1080-1004);
    left_click(Mouse);
end % end function

function start_capture(Mouse)
    % move mouse to record button (when stem fullscreen)
    Mouse.mouse.mouseMove(1720, 1080-170);
    left_click(Mouse);
end % end function

function stop_capture(Mouse, isCorrect, response)
    % called with isCorrect= 1 if stimulus should have elicitied user
    % response
    
    % move mouse to record button (when stem fullscreen)
    Mouse.mouse.mouseMove(1720, 1080-170);
    left_click(Mouse);
    
    % save file
    name= horzcat('was_', isCorrect, '_ans_', response);
    clipboard('copy', name);
    clipboard('paste');
    % click confirm button
    Mouse.mouse.mouseMove(1720, 1080-493);
    left_click(Mouse);
end % end function

function left_click(Mouse)    
    Mouse.mouse.mousePress(Mouse.left_click);
    Mouse.mouse.mouseRelease(Mouse.left_click);
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
