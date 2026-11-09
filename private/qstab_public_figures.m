function C=qstab_public_figures(F,root)
%QSTAB_PUBLIC_FIGURES Copy this run's exported PNG/PDF bytes into the MAT report.
% Does not render, recompute, modify or delete images. Validates origin, export
% metadata and signatures; the resulting files still await visual inspection.
% Only <root>/qstab_paper_figures/<id>.png/.pdf is accepted. No arbitrary file
% attachments. A conservative 64 MiB total cap keeps the existing MAT-v7 report.
if ~isstruct(F)||~isscalar(F)||~all(isfield(F,{'ok','items','export'}))||~F.ok|| ...
        ~isfield(F.export,'ok')||~F.export.ok||~isfield(F.export,'files')||isempty(F.items)
    error('qstab:PublicFigures','Completed figure data and export are required.');
end
C=struct('schema','qstab.figure-review.1','status','CAPTURED','ok',true, ...
    'visualInspectionPerformed',false,'totalBytes',0, ...
    'files',struct('id',{},'format',{},'name',{},'bytes',{},'sha256',{},'data',{}), ...
    'note','Exact MATLAB export bytes, not Python redraws; signatures do not establish visual quality.');
ids={F.items.id};formats={'png','pdf'};
if numel(unique(ids))~=numel(ids),error('qstab:PublicFigures','Duplicate figure IDs.');end
for j=1:numel(ids)
    id=ids{j};
    if ~ischar(id)||isempty(regexp(id,'^[A-Za-z0-9_-]+$','once'))
        error('qstab:PublicFigures','Unsafe figure ID.');
    end
    for k=1:numel(formats)
        fmt=formats{k};expected=fullfile(root,'qstab_paper_figures',[id '.' fmt]);
        ix=find(strcmp({F.export.files.id},id)&strcmp({F.export.files.format},fmt));
        if numel(ix)~=1||~localSame(F.export.files(ix).file,expected)
            error('qstab:PublicFigures','Missing, duplicated or out-of-directory %s.%s.',id,fmt);
        end
        d=dir(expected);
        if numel(d)~=1||d.isdir||d.bytes~=F.export.files(ix).bytes||d.bytes<8|| ...
                C.totalBytes+d.bytes>64*1024^2
            error('qstab:PublicFigures','Missing, changed, empty or excessive export: %s.',expected);
        end
        bytes=localRead(expected);
        if numel(bytes)~=d.bytes,error('qstab:PublicFigures','File size changed during read.');end
        if strcmp(fmt,'png'),sig=uint8([137 80 78 71 13 10 26 10]);
        else,sig=uint8('%PDF-');end
        if ~isequal(bytes(1:numel(sig)).',sig),error('qstab:PublicFigures','Wrong export signature.');end
        digest='';
        try
            md=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
            md.update(typecast(bytes(:),'int8'));h=typecast(md.digest(),'uint8');
            digest=lower(reshape(dec2hex(h,2).',1,[]));
        catch
            % Empty digest explicitly records that no digest was obtained.
        end
        C.files(end+1)=struct('id',id,'format',fmt,'name',[id '.' fmt], ...
            'bytes',numel(bytes),'sha256',digest,'data',bytes); %#ok<AGROW>
        C.totalBytes=C.totalBytes+numel(bytes);
    end
end
end
function b=localRead(p)
fid=fopen(p,'rb');if fid<0,error('qstab:PublicFigures','Cannot read export.');end
closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
b=fread(fid,Inf,'*uint8');
end
function tf=localSame(a,b)
a=strrep(a,'\','/');b=strrep(b,'\','/');
if ispc,tf=strcmpi(a,b);else,tf=strcmp(a,b);end
end
