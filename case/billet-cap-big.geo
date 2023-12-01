SetFactory("OpenCASCADE");

xc = 0.0;
yc = 0.78;
zc = 0.0;
outerboxh = 0.6; obh = outerboxh/2; //outerboxh = 0.2;  obh = outerboxh/2;
innerboxh = 0.6; ibh = innerboxh/2; //innerboxh = 0.08; ibh = innerboxh/2;
rad = 1;

Box(1) = {xc-2*obh,yc-0.5,zc-2*obh, 2*outerboxh, 2.2*outerboxh, 2*outerboxh};
Box(2) = {xc-  ibh,yc-0.1,zc-  ibh,   innerboxh,        0.62,   innerboxh};
Sphere(3) = {0.0, 0.0, 0.0, rad};
V1 = BooleanDifference { Volume{1}; Delete;} { Volume{2:3}; Delete;};

Transfinite Line {25:36} = 31 Using Progression 1.0;
Transfinite Line {15,19,23,24,37:44} = 21 Using Progression 1.0;

Transfinite Surface {10:14,16:20};
Recombine Surface {10:20};
Recombine Volume {V1};

Physical Volume ("Volume") = {V1};
Physical Surface ("INNER-BOX") = {10,17:20};
Physical Surface ("SPHERE") = {15};
Physical Surface ("OUTER-BOX") = {11:16};
