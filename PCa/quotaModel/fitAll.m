function fitAll(fixdata)

if nargin==0
    fixdata = 0;
end

parfor j = 2:7
    fitCase(j,fixdata);
end

seeParams;
for j=1:7
    plotPsa(j,0,1,0);
end

end
