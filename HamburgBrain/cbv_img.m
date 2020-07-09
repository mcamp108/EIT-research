function imgOut = cbv_img(pig)
% -------------------------------------------------------------------------
% DESCRIPTION:
%
% somehow need to find where in MRI image is the brain
% then need to know where brain is in EIT images
% then need to compare center of gravity (foci location) in each
% then need to compare them to see how close / far they are.
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

imgDir = 'C:\Users\Mark\Documents\GraduateStudies\LAB\HamburgBrain\Figures\paper\cbv\';

switch pig
    case '8-2'   
        imgName     = '8-2_Perfusion_Weighted_split_1_z32.png';
        brainTop    = 30;
        brainBot    = 365;
        brainLeft   = 525;
        brainRight  = 860;
    case '9-2'     
        imgName     = '9-2_SWI - as a 30 frames MultiVolume by AcquisitionTime_z61.png';
        brainTop    = 25;
        brainBot    = 350;
        brainLeft   = 501;
        brainRight  = 835;
    case '10-2'    
        imgName     = '10-2_RELCBV_LOCAL_ep2d_perf_p2_MoCo_1_z70.png';
        brainTop    = 20;
        brainBot    = 360;
        brainLeft   = 510;
        brainRight  = 865;
    case '11-2'    
        imgName     = '11-2_RELCBV_LOCAL_ep2d_perf_p2_MoCo_z67.png';
        brainTop    = 35;
        brainBot    = 400;
        brainLeft   = 490;
        brainRight  = 870;        
    case '12-2'    
        imgName     = '12-2_RELCBV_LOCAL_ep2d_perf_p2_MoCo_z57.png';
        brainTop    = 55;
        brainBot    = 345;
        brainLeft   = 495;
        brainRight  = 825;
end % end switch

img1 = imread(horzcat(imgDir, imgName));
imgOut = double( img1(brainTop:brainBot,brainLeft:brainRight) );

end % end function