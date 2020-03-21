function name = caller()

stack = dbstack;

name = '';
if length(stack) > 1
   caller = stack(2).file;
end