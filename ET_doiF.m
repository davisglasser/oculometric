function [directions results raw_position RTs] = ET_doiF(doi, durations, raw_position, results, directions, RTs)

count=0;
dir_doi = [];
result_doi = dir_doi;
pos_doi = dir_doi;
for i=1:length(directions)
    if round(durations(i))==round(doi)
        count = count+1;
        dur_doi(count)=durations(i);
        pos_doi(count,:)=raw_position(i,:);
        result_doi(count)=results(i);
        dir_doi(count)=directions(i);
        RTs_doi(count)=RTs(i);
    end
end
directions = dir_doi;
results = result_doi;
raw_position = pos_doi;
RTs = RTs_doi;
samp = 2;