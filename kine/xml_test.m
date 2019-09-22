%xmlµÄ¶ÁÈ¡
% D1=parseXML('input.xml');
% D2=xmlread('Epson_G10.xml','r');
p560xml = parseXML('puma560.xml');
% D2=parseXML('Epson_G10.xml');
type=p560xml.Name;
name=p560xml.Attributes;
if strcmp(name.Name,'name')
    temp.name=name.Value;
end
children=p560xml.Children;
n=length(children);
k=0;
for i=1:n
    if ~isempty(children(i).Attributes)
        k=k+1;
        child2(k)=children(i);
    end
end
n_child=length(child2);
n_robot=n_child-3;
base=child2(1);
tool=child2(n_child-1);
T0=child2(n_child);

