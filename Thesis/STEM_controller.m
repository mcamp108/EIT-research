function STEM_controller(varargin)
% -------------------------------------------------------------------------
% DESCRIPTION:
%   STEM_controller(varargin) 
%
%   Use this function to control the Swisstom EIT
%   Monitor. This function is capable of launching STEM, starting a STEM
%   recording, stopping (and saving) a STEM recording, or executing a left
%   click on a java robot mouse.
%
% -------------------------------------------------------------------------
% PARAMETERS:
%   Mouse:
%       
%   command:
%       one of "launch stem", "left click", "start capture", or "stop
%       capture"
%   example usages:
%       STEM_controller(Mouse, "launch stem")
%       STEM_controller(Mouse, "left click")
%       STEM_controller(Mouse, "start capture")
%       STEM_controller(Mouse, "stop capture", stimulus_time, response_time, answers, responses)
% -------------------------------------------------------------------------   
% RETURNS:
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

if command== "launch stem"
    launch_stem(Mouse);
elseif command== "left click"
    left_click(Mouse);
elseif command== "start capture"
    start_capture(Mouse);
elseif command== "stop capture"
    stimulus_time= varargin{3};
    response_time= varargin{4};
    answers= varargin{5};
    responses= varargin{6};
    stop_capture(Mouse, stimulus_time, response_time, answers, responses);
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


function stop_capture(Mouse, stimulus_time, response_time, answers, responses)

% move mouse to record button (when stem fullscreen)
Mouse.mouse.mouseMove(1720, 1080-170);
left_click(Mouse);

% save file
name= char(date);
clipboard('copy', name);
clipboard('paste');
% click confirm button
Mouse.mouse.mouseMove(1720, 1080-493);
left_click(Mouse);

% save timing data
header= {'stimulus time', 'response time', 'answers', 'responses'};
A= [stimulus_time, response_time, answers, responses];
save_name= horzcat(name, '.csv');
xlswrite(save_name, header, 'A1' );
xlswrite(save_name, A, 'A2' );

end % end function

