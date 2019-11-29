
// Trilateration example
// from Wikipedia
// 
// pA, pB and pC are the centres of the spheres
// If necessary the spheres must be translated
// and rotated so that:
// -- all z values are 0
// -- pA is at the origin
pA = [0,0,0];
// -- pB is on the x axis
pB = [10,0,0];
pC = [7,7,0];

// rA , rB and rC are the radii of the spheres
rA = 8;
rB = 5;
rC = 7;

if ( pA != [0,0,0]){
   echo ("ERROR: pA must be at the origin");
   assert(false);
}

if ( (pB[2] !=0 ) || (pC[2] !=0)){
   echo("ERROR: all sphere centers must be in z = 0 plane");
   assert(false);
}

if (pB[1] != 0){
   echo("pB centre must be on the x axis");
   assert(false);
}

// show the spheres
module spheres(){
   translate (pA){
      sphere(r= rA, $fn = rA * 10);
   }

   translate(pB){
      sphere(r = rB, $fn = rB * 10);
   }

   translate(pC){
      sphere (r = rC, $fn = rC * 10);
   }
}

function unit_vector( v) = v / norm(v);

ex = unit_vector(pB - pA) ;
echo(ex = ex);

i = ex * ( pC - pA);
echo (i = i);

ey = unit_vector(pC - pA - i * ex);
echo (ey = ey);

d = norm(pB - pA);
echo (d = d);

j =  ey * ( pC - pA);
echo (j = j);

x = (pow(rA,2) - pow(rB,2) + pow(d,2)) / (2 * d);
echo( x = x);

// size of the cube to subtract to show 
// the intersection of the spheres
cube_size = [10,10,10];

if ( ((d - rA) >= rB) || ( rB >= ( d + rA)) ){
   echo ("Error Y not solvable");
   spheres();
}else{
   y = (( pow(rA,2) - pow(rC,2) + pow(i,2) + pow(j,2)) / (2 * j))
      - ( i / j) * x;
   echo(y = y);
   zpow2 = pow(rA,2) - pow(x,2) - pow(y,2);
   if ( zpow2 < 0){
      echo ("z not solvable");
      spheres();
   }else{
      z = sqrt(zpow2);
      echo (z = z);
      // subtract a cube with one of its corners 
      // at the point where the sphers intersect
      difference(){
         spheres();
         translate ([x,y - cube_size[1],z]){
           cube(cube_size);
         }
      }
      translate ([x,y - cube_size[1],z]){
           %cube(cube_size);
      }
  }
}
