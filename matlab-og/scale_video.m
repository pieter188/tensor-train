%% mask function for frames of video
% Pieter van Klaveren
% 24-01-2021
% input: Video via VideoReader, desired resolution
% output: scaled video (uint8)

function scaled_Video = scale_video(Video,res)
% depending on the type of input -> convert to uint8
if strcmp(class(Video),'VideoReader') == 1
    % read all frames -> RBG data of each frame in 4D tensor
    frames = read(Video,[1 Inf]);
elseif strcmp(class(Video),'uint8') == 1
        frames = Video;
end
% obtain size of video
a = size(frames);

for i=1:a(4)
    scaled_Video(:,:,:,i) = imresize(frames(:,:,:,i),res);
end

end

