

#include <quan/out/length.hpp>
#include <quan/two_d/out/vect.hpp>
#include <quan/fixed_quantity/literal.hpp>
#include <cassert>

/*
  Trilateration example from Wikipedia https://en.wikipedia.org/wiki/Trilateration

  also see the Openscad example which really helps to visualise the problem in 3D
  https://github.com/kwikius/Trilateration/blob/master/trilateration.scad

  https://stackoverflow.com/questions/9747227/2d-trilateration
  https://stackoverflow.com/questions/23400351/localizing-a-point-using-distances-to-three-other-points-in-3-d/23401529

  requires my quan library ( headers only required)
  https://github.com/kwikius/quan-trunk
*/

namespace {
    QUAN_QUANTITY_LITERAL(length,km)
    typedef quan::two_d::vect<quan::length::km> point;
    typedef quan::two_d::vect<quan::length::km> km_vect;
    typedef quan::two_d::vect<double> dvect;
}


int main()
{

   // we use 2d points here. 3D points could be used, but their z values must be 0
   // point pA must be at the origin
   // TODO
   // Translate by V all points so that pA is at origin
   // rotate  all points by alpha until pB is on the x axis
   // do calc
   // rotate result by - alpha
   // translate result by -V
   point constexpr pA{0.0_km,0.0_km};

   assert( (pA == point{0.0_km, 0.0_km}));

   // point pB must be on the x axis
   point constexpr pB{10.0_km,0.0_km};

   assert( pB.y == 0.0_km);

   point constexpr pC{7.0_km,7.0_km};

   auto constexpr rA = 8.0_km;
   auto constexpr rB = 5.0_km;
   auto constexpr rC = 7.0_km;

   auto const ex = unit_vector(pB-pA);

   auto const i = dot_product(ex,(pC- pA));
   
   auto const ey = unit_vector(pC - pA - i * ex) ;

   auto const d = magnitude(pB - pA);

   auto const j = dot_product(ey,pC - pA);
    
   std::cout << "ex = " << ex << '\n';
   std::cout << "i  = " << i  << '\n';
   std::cout << "ey = " << ey << '\n';
   std::cout << "d  = " << d  << '\n';
   std::cout << "j  = " << j  << '\n';

   auto const x = (quan::pow<2>(rA) - quan::pow<2>(rB) + quan::pow<2>(d)) / ( 2 * d);
   std::cout << "x = " << x << '\n';

   if ( ((d - rA) >= rB ) || (rB >= (d + rA))){
      std::cout << "y : no solution\n";
      return 0;
   }

   auto const y =  (
         ( quan::pow<2>(rA) - quan::pow<2>(rC) + quan::pow<2>(i) + quan::pow<2>(j))
               / ( 2 * j) 
                  ) - ( i / j) * x;
   
   std::cout << "y = " << y << '\n';

   auto z_2 = quan::pow<2>(rA) - quan::pow<2>(x) - quan::pow<2>(y);
   if ( z_2 >= quan::pow<2>(0_km)){
      auto const z = sqrt(z_2);
      std::cout << "z = " << z << '\n';
   }else{
      std::cout << "z : no solution\n";
   }
}
