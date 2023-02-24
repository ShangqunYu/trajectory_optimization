function D = fcn_get_D(in)

d = in(1);
e = in(2);
f = in(3);
D = [0 0 0;
   e -d 0;
   f 0 -d;
   -e d 0;
   0 0 0;
   0 f -e;
   -f 0 d;
   0 -f e;
   0 0 0];
end