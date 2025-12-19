function varargout = pickspot(inputimage)
% varargout = pickspot(inputimage)
 
if length(size(inputimage))>3
  error('only sizes up to 3 dimensions are supported')
end

if length(size(inputimage))==3
  x=squeeze(inputimage(:,:,0));
else
  x=squeeze(inputimage);
end

coords = hand_pick(x);
co=newim(x);
diptruesize(30);
for jj=1:size(coords,1)
    co(coords(jj,1),coords(jj,2))=1;
end
co=overlay(stretch(x),co==1);
varargout{1}=coords;
varargout{2}=co;