%% rgb2gray for videos function
% Pieter van Klaveren
% 24-01-2021
% input: Video in rgb format, animation on or off
% output Video in grayscale format (uint8)

function gray_Video = rgb2gray_video(Video)
% depending on the type of input -> convert to uint8
if strcmp(class(Video),'VideoReader') == 1
    % read all frames -> RBG data of each frame in 4D tensor
    frames = read(Video,[1 Inf]);
elseif strcmp(class(Video),'uint8') == 1
        frames = Video;
end

% obtain size of the frames
a = size(frames);
% reshape so frames are aligned with vertical pixels (row of tensor)
frames = permute(frames,[1,4,2,3]);
% combine amount of vertical pixels and frames
frames = reshape(frames,[a(1)*a(4),a(2),a(3)]);
% apply colour to grayscale conversion
frames = rgb2gray(frames); % must have uint8 type as input
% substracting frames from rows
frames = reshape(frames,[a(1),a(4),a(2),1]);
% right order of rows, column, gray, frames
gray_Video = permute(frames,[1,3,4,2]);
% fprintf('Video was in RBG format, but is transformed to grayscale format.\n')

end