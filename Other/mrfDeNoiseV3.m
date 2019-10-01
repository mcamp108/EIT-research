function [imageDenoised, betheE]= mrfDeNoiseV3(image, lambda)
% This function was written by :
%                               Mark Campbell
%                               Carleton University
imageRnd= round(image);
numLabels= max(imageRnd(:))- min(imageRnd(:))+ 1; % Number of possible intensity values
[m, n]= size(image);
pixels= m* n; % Number of pixels in input image
nonNaN= find( ~isnan( image))';
meanImage= mean(image(nonNaN));

A = getAdjMat(image); % First build adjacency matrix where value of edge is edge index
[i,j,e] = find(triu(A)); % i is row, j is column, e is edge number
numEdges = numel(e); % Number of unique edges. Messages will be passed in forward and reverse directions acoss the same edge

nodeBeliefs= zeros(numLabels, pixels); % Stores variable node potentials. numLabels x pixels
nodePotentials= zeros(numLabels, pixels); % Contains all possible unary potential values
edgePotentials= zeros(numLabels, numLabels); % Contains all possible binary potential values
betheE= []; % Vector to store Bethe energy output for each iteration

factNodesPrev= zeros(4, numLabels, pixels); % Stores a nodes' messages from its neighbours for time = t-1
factNodes= zeros(4, numLabels, pixels); % Stores a nodes' messages from its neighbours for time = t. First index is direction of message relative to recipient ndoes [up, right, down, left]
validEdges= zeros(4, pixels); % Store valid edges in case messages become 0

% image to be used in calculations (integer values, minimum value brought to 1
imageCopy= imageRnd(nonNaN) + abs(min(image(:)))+ 1;

iter= 0;
doneIter= false;
done= false;
cycle= 1;
while done== false
    iter= iter+ 1;
    disp("Cycle: "+ cycle+ " Iteration: " + iter);
    
    for row= 1:numLabels
%         edgePotentials(row, :)= exp(- (lambda.* min(abs(row- (1:numLabels)), 0.75)) ); % truncated quadtratic smoothness cost function
%         edgePotentials(row, :)= (1+ 1/2*( (row- (1:numLabels) )/ sigma).^2).^ (-lambda); % based on student's t-distribution (Lan et al.)
        edgePotentials(row, :)= exp(- (abs((row- (1:numLabels))).* lambda)); % standard smoothness cost function
        nodePotentials(row, nonNaN)= exp(- ((row- imageCopy(:)).^2) );
    end % end for

% 1. Variable node sends message to factor node that is the product of the messages recieved from its neighbouring factor nodes other than receiver
for edge= 1: numEdges % for each edge
    sender= i(edge); % pixel number of message sender
    receiver= j(edge); % pixel number of message recipient
    [fwd, rev] = messageDirection(sender, receiver, m); % factNode index of receiver corresponding to cardinal position of sender relative to receiver (see messageDirection function below)
    
    if iter== 1
       factNodes(fwd, :, receiver)= 1; % Initialize receiver factor potential to 1 in direction fwd for all labels (uninformed prior).
       factNodes(rev, :, sender)= 1; % Initialize sender factor potential to 1 for all labels in direction rev
       validEdges(fwd, receiver)= 1; % Has edge with node= 1, does not have edge= 0
       validEdges(rev, sender)= 1;
    else
        for dir = [fwd, rev]
            revdir= mod(dir+2, 4);
            validEdgeIndex= find(validEdges(:, sender)); % Indexes of non-zero entries specify direction valid edges
            validEdgeIndex= validEdgeIndex(validEdgeIndex~= revdir); % Excludes message received from receiver
            localEvidence= edgePotentials.* prod(factNodesPrev(validEdgeIndex, :, sender), 1)'; % Multiply along rows
            message= sum(nodePotentials(:, sender).* localEvidence, 1);
            factNodes(dir, :, receiver)= message;
            temp= sender; % Prepare to pass message in the reverse direction
            sender= receiver;
            receiver= temp;
        end % end for
    end % end if
end % end for

% 2. Normalize messages received by a neighbouring node across all labels
factNodes= normalize_messages(factNodes, nonNaN); 
factNodesPrev= factNodes; % Move to the next iteration. Messages(t) -> Messages(t-1) relative to next time frame

% Bethe free energy calculation components
Exy= 0; % First term.
Ex= 0; % Second term.
Hx= 0; % Third term.
Ixy= 0; % Fourth term.
nodeBeliefPrev= nodeBeliefs;

% 3. Calculate node beliefs.
nodeBeliefPrev= nodeBeliefs;
[nodeBeliefs, Ex, Hx]= calc_nodeBeliefs(nodeBeliefs, nodePotentials, validEdges, factNodes, nonNaN);

% 4. Calculate edge beliefs.
[Exy, Ixy]= calc_edgeBeliefs(i, j, m, numEdges,  numLabels, validEdges, factNodes, nodePotentials, edgePotentials);

% 5. Calculate Bethe free energy (Yedidia 2001) and terminate function if at a minimum
if iter> 1
    betheFreeE= -Exy+ Ex+ Ixy- Hx;
    if iter>2 && (betheFreeEPrev < betheFreeE || betheFreeEPrev/ betheFreeE < 1.05 || betheFreeEPrev== 0)
        nodeBeliefs= nodeBeliefPrev;
        doneIter= true;
        betheE= [betheE; betheFreeEPrev];
    else
        betheFreeEPrev= betheFreeE;
    end % end if
end % end if

% 6. Select most probable label for each pixel once beliefs converge
if doneIter 
    if cycle== 1 || (cycle>= 2 && (betheE(end) > betheE( length(betheE)-1) ))
        cycle= cycle+ 1;
        iter= 0;
        doneIter= false;
        imageDenoised= get_map_estimate(image, nonNaN, nodeBeliefs);
        imageCopy= imageDenoised(nonNaN);
    else % If this loop did not produce a lower betheE, return the image from the previous iteration
        done= true;
        betheE= betheE(1: length(betheE)-1);
    end % end if
end % end if
end % end while (iter)

meanImageDenoised= mean(imageDenoised(nonNaN));
imageDenoised= imageDenoised- (meanImageDenoised- meanImage); % keep mean of output image same as mean of input image

end % end function


function [dir, revdir] = messageDirection(i, j, numRows)
% parameters: sender int i, recipient int j, and int number of rows in the image.
% returns: direction in which message is being received by recipient 
% 1= north, 2= east, 3= south, 4= west
result= j-i;
if result== 1
    dir= 1; % north
elseif result == -numRows
    dir= 2; % east
elseif result == -1
    dir= 3; % south
else
    dir= 4; %west
end % end if
revdir = mod(dir+2, 4); % reverse message
end % end function

function factNodes= normalize_messages(factNodes, nonNaN)
for pix= nonNaN
    for dir= 1:4
        Messages= factNodes(dir, :, pix);
        Z= sum(Messages);
        if Z~= 0
            factNodes(dir, :, pix) = Messages/Z;
        end % end if
    end % end for
end % end for

end % end function

function [nodeBeliefs, Ex, Hx]= calc_nodeBeliefs(nodeBeliefs, nodePotentials, validEdges, factNodes, nonNaN)
Ex= 0;
Hx= 0;
for pixel= nonNaN % update belief for each pixel
    validEdgeIndex= find(validEdges(:, pixel));
    nodeBeliefs(:, pixel)= nodePotentials(:, pixel).* prod(factNodes(validEdgeIndex, :, pixel), 1)'; % Overall probability of pixel having each label
    nonZeroNodeBeliefIndex= find(nodeBeliefs(:, pixel));
    nonZeroNodeBeliefs= nonzeros(nodeBeliefs(:, pixel));
    Q= length(validEdgeIndex);
    Ex= Ex+ (Q* dot(nonZeroNodeBeliefs, log(nodePotentials(nonZeroNodeBeliefIndex, pixel))));
    Hx= Hx+ (Q* dot(nonZeroNodeBeliefs, log(nonZeroNodeBeliefs)));
end % end for
end % end function

function [Exy, Ixy]= calc_edgeBeliefs(i, j, m, numEdges,  numLabels, validEdges, factNodes, nodePotentials, edgePotentials)
Exy= 0;
Ixy= 0;
edgeBeliefs= zeros(numLabels, numLabels, numEdges); % Stores probability associated with each label permutation for each pair of pixels
for edge= 1: numEdges
    Xi= i(edge); % Pixel number of message sender
    Xj= j(edge); % Pixel number of message recipient
    [fwd, rev] = messageDirection(Xi, Xj, m); % factNode index of receiver corresponding to cardinal position of sender relative to receiver (see messageDirection function below)
    validEdgeIndexFwd= find(validEdges(:, Xi));
    validEdgeIndexFwd= validEdgeIndexFwd(validEdgeIndexFwd~= rev); % Excludes message received from receiver
    XiIncomingMessages= prod(factNodes(validEdgeIndexFwd, :, Xi), 1)';
 
    validEdgeIndexRev= find(validEdges(:, Xj));
    validEdgeIndexRev= validEdgeIndexRev(validEdgeIndexRev~= fwd); % Excludes message received from receiver
    XjIncomingMessages= prod(factNodes(validEdgeIndexRev, :, Xj), 1)';
    
    localPotentials= nodePotentials(:,Xi).* nodePotentials(:,Xj).* edgePotentials;

    edgeBeliefs(:, :, edge)= localPotentials.* XiIncomingMessages.* XjIncomingMessages;
    nonZeroEdgeBeliefIndex= find(edgeBeliefs(:, :, edge));
    nonZeroEdgeBeliefs= nonzeros(edgeBeliefs(:, :, edge));
    Exy= Exy+ dot( nonZeroEdgeBeliefs, log(localPotentials(nonZeroEdgeBeliefIndex)) );
    Ixy= Ixy+ dot( nonZeroEdgeBeliefs, log(nonZeroEdgeBeliefs) );
end % end for
end % end function

function imageDenoised= get_map_estimate(image, nonNaN, nodeBeliefs)
imageDenoised= image;
for pixel= nonNaN
    label= find( nodeBeliefs(:, pixel)== max(nodeBeliefs(:, pixel)) );
    if length(label) > 1
        label= label(randi([1,length(label)])); % Choose randomly if there is a tie.
    end
    imageDenoised(pixel)= label; % Choose the label with the highest probability
end % end for
end % end function