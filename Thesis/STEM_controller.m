function STEM_controller(varargin)

% -------------------------------------------------------------------------
% DESCRIPTION:
%
%
%
% -------------------------------------------------------------------------
% PARAMETERS:
% 
%
%
% -------------------------------------------------------------------------   
% RETURNS:
% 
%
%
% -------------------------------------------------------------------------   
% AUTHOR:
%   Mark Campbell
%   Carleton University
%   markacampbell@cmail.carleton.ca
% -------------------------------------------------------------------------
% VERSION:
%   1.0.0
% -------------------------------------------------------------------------
Mouse= varargin{1};
command= varargin{2};

% record mouse's original location
C = get(0, 'PointerLocation');

if nargin== 4
    isCorrect= varargin(3);
    response= varargin(4);
end % end if

if command== "launch stem"
    launch_stem(Mouse);
elseif command== "left click"
    left_click(Mouse);
elseif command== "start capture"
    start_capture(Mouse);
elseif command== "stop capture"
    if exist(isCorrect, 'var') && exist(response, 'var')
        stop_capture(Mouse, isCorrect, response);
    else
        disp("Missing arguments argin.isCorrect and/or argin.response.");
    end % end if
end

% move mouse back to its original location
Mouse.mouse.mouseMove(C(1), 1080-C(2));

end % end function

function launch_stem(Mouse)

cd 'C:\Users\Mark\Desktop\STEM';
system('run.bat');
cd 'C:\Users\Mark\Documents\GraduateStudies\LAB\Thesis\Data\Headcap';
mkdir(date);
% click sensor quality screen
Mouse.mouse.mouseMove(232, 1080-1004);
left_click(Mouse);

end % end function

function left_click(Mouse)    
 
Mouse.mouse.mousePress(Mouse.left_click);
Mouse.mouse.mouseRelease(Mouse.left_click);
    
end % end function

function start_capture(Mouse)

% move mouse to record button (when stem fullscreen)
Mouse.mouse.mouseMove(1720, 1080-170);
left_click(Mouse);
    
end % end function

function stop_capture(Mouse, isCorrect, response)

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

