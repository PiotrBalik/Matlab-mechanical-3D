warning('non-square basis is not tested (e.g. 3x4, 4x6...')
x=3;
y=3;

storey=6;

%per storey:
ncol = (x+1)*(y+1);
nbeam = (x+1)*y+x*(y+1);
NS=(ncol+nbeam);
nele = NS*storey;
njoint = ncol*storey;
NEL=[1:nele];
type = mod(NEL-1,ncol+nbeam)+1;
type(type<=ncol)=2;
type(type~=2)=1;


%element No.   node from/to  sec. type
eltable = [ NEL' zeros(nele,2) type'];

fromC=1:2:ncol*2;
toC=2:2:ncol*2;
connCol=[fromC' toC'];

connBeam=[toC(1:end-1)' toC(2:end)'];
idxs=(x+1):(x+1):(x+1)*y;
connBeam(idxs,:)=[];

fromB=2:(x+1)*(y-1):ncol*2;
fromB=fromB'+[0:2:(x+1)*(y-1)];

from=fromB(1:end-1,1:x+1);
to=fromB(2:end,1:x+1);
connAdd=[from(:) to(:)];

connBeam=[connBeam; connAdd];

eltable(1:ncol,2:3)=connCol;
eltable(ncol+1:NS,2:3)=connBeam;

for id=2:storey
    
    fromC=toC;
    toC=(ncol*id+1):ncol*(id+1);
    connCol=[fromC' toC'];

    connBeam=[toC(1:end-(y+1))' toC(y+2:end)'];

    connAdd=[toC(1:end-1)' toC(2:end)'];
    connAdd(idxs,:)=[];

    connBeam=[connBeam; connAdd]; %#ok

    eltable( (1+NS*(id-1)):(NS*(id-1)+ncol),2:3)=connCol;
    eltable( (NS*(id-1)+ncol+1):NS*id,2:3)=connBeam;

end


fromC=1:2:ncol*2;
toC=2:2:ncol*2;
%group assignment
A = [fromC(1:2:2*(x+1)) toC(2:2:2*(x+1))];
vec = [0:storey-1]'*40;
A = [A+vec];
B = A(storey/2+1:end,:);
B = sort(B(:))';
A = A(1:storey/2,:);
A=sort(A(:))';

C = [fromC(2:2:2*(x+1)) toC(1:2:2*(x+1))];
vec = [0:storey-1]'*40;
C = [C+vec];
D = C(storey/2+1:end,:);
D = sort(D(:))';
C = C(1:storey/2,:);
C=sort(C(:))';

F = [(nele-nbeam+1):nele];

E = setdiff(NEL,[A B C D F]);
letters=[ismember(NEL,A)*'A'+ismember(NEL,B)*'B'+ismember(NEL,C)*'C'+ismember(NEL,D)*'D'+...
    ismember(NEL,E)*'E'+ismember(NEL,F)*'F'];
fprintf('\nEl no.\tfrom\tto\t\ttype\tgroup\n')
fprintf([repmat('%d\t\t',[1 4]) '%c\n'],[eltable letters' ]')
