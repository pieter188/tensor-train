%% mask function for frames of video
% Pieter van Klaveren
% 24-01-2021
% input: video via VideoReader, percentage of corrupted pixels values between
%        0 and 1.
% output: masked video (uint8) with a fraction alpha of the pixels missing
function masked_Video = mask_video(Video,alpha,beta)
%% obtain video information
% depending on the type of input -> convert to uint8
if strcmp(class(Video),'VideoReader') == 1
    % read all frames -> RBG data of each frame in 4D tensor
    frames = read(Video,[1 Inf]);
elseif strcmp(class(Video),'uint8') == 1
        frames = Video;
end

% from unit8 to double
frames = double(frames);

% take away initial frames that are not masked
a = size(frames);
masked_Video = zeros(a);
masked_Video(:,:,:,1:round(beta*a(4))) = frames(:,:,:,1:round(beta*a(4)));
% creating mask
mask = rand([size(frames(:,:,:,round(a(4)*beta):end))]);
mask(mask>alpha) = 1;
mask(mask<=alpha) = 0;
masked_Video(:,:,:,round(a(4)*beta):end) = mask.*frames(:,:,:,round(a(4)*beta):end);
% from double to uint8
masked_Video = uint8(masked_Video);
end