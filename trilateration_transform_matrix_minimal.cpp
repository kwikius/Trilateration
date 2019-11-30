
#include <cstdlib>
#include <iostream>

#include <quan/out/angle.hpp>
#include <quan/atan2.hpp>
#include <quan/out/length.hpp>
#include <quan/three_d/out/vect.hpp>
#include <quan/three_d/rotation.hpp>
#include <quan/three_d/sphere.hpp>

#include <quan/fusion/make_matrix.hpp>
#include <quan/fusion/matrix.hpp>
#include <quan/fusion/make_translation_matrix.hpp>
#include <quan/fusion/make_3d_x_rotation_matrix.hpp>
#include <quan/fusion/make_3d_y_rotation_matrix.hpp>
#include <quan/fusion/make_3d_z_rotation_matrix.hpp>
#include <quan/fusion/make_column_matrix.hpp>
#include <quan/fun/display_matrix.hpp>
#include <fstream>

// calc diagnostic output
#define DEBUG_PRINT

#define SHOW_MATRIX_CALC

namespace {

   QUAN_QUANTITY_LITERAL(length,km)
   typedef quan::three_d::vect<quan::length::km > point;

   typedef quan::three_d::sphere<quan::length::km> sphere;

#if defined DEBUG_PRINT
   std::ostream & operator<< ( std::ostream & out, sphere const & c)
   {
      return out << "sphere(centre = " << c.centre << ", radius = " << c.radius << ")";
   }
#endif

   auto constexpr epsilon_km = 1.e-6_km;

// to remove unnecessary calcs ( otherwise useful for exposition )
#define MINIMAL_CALCS
   // A B C must be normalised
   // where A is centred at origin
   // B is centred on x axis
   // C is centred on xy plane
   bool ll_trilaterate( sphere const& A, sphere const & B, sphere const & C, point & intersection_point)
   {
      assert( (A.centre == point{0.0_km, 0.0_km,0.0_km}));
      assert(abs(B.centre.y) < epsilon_km); 
      assert(abs(B.centre.z) < epsilon_km);
      assert(abs(C.centre.z) < epsilon_km); 

#if defined MINIMAL_CALCS
      auto const ex = unit_vector(B.centre);   // direction of B to origin
      auto const i = dot_product(ex,(C.centre));   
      auto const ey = unit_vector(C.centre - i * ex) ;
      auto const d = magnitude(B.centre);       // distance B to origin
      auto const j = dot_product(ey,C.centre);
#else
      auto const ex = unit_vector(B.centre-A.centre);
      auto const i = dot_product(ex,(C.centre- A.centre));
      auto const ey = unit_vector(C.centre - A.centre - i * ex) ;
      auto const d = magnitude(B.centre - A.centre);
      auto const j = dot_product(ey,C.centre - A.centre);
#endif
      auto const x = (quan::pow<2>(A.radius) - quan::pow<2>(B.radius) + quan::pow<2>(d)) / ( 2 * d);

      if ( ((d - A.radius) >= B.radius ) || (B.radius >= (d + A.radius))){
         std::cout << "y : no solution\n";
         return false;
      }

      auto const y =  (
            ( quan::pow<2>(A.radius) - quan::pow<2>(C.radius) + quan::pow<2>(i) + quan::pow<2>(j))
                  / ( 2 * j) 
                     ) - ( i / j) * x;

      auto const z_2 = quan::pow<2>(A.radius) - quan::pow<2>(x) - quan::pow<2>(y);
      if ( z_2 >= quan::pow<2>(0_km)){
         auto const z = sqrt(z_2);

         intersection_point = point{x,y,z};
   #if defined DEBUG_PRINT
         std::cout << "ll intersection point = " << intersection_point << '\n';
   #endif
         return true;
        
      }else{
         std::cout << "z : no solution\n";
         return false;
      }
      return true;
   }

}

bool trilaterate(sphere const& A, sphere const & B, sphere const & C,point & out)
{
   auto const pA0 = A.centre;
   auto const pB0 = B.centre;
#if defined SHOW_MATRIX_CALC
   auto  pB0v = quan::fusion::make_column_matrix(B.centre);
#endif

   auto const pC0 = C.centre;
#if defined DEBUG_PRINT
   std::cout << "pA0 = " << pA0 << '\n';
   std::cout << "pB0 = " << pB0 << '\n';
#if defined SHOW_MATRIX_CALC
   display(pB0v,"pB0v = ");
#endif
   std::cout << "pC0 = " << pC0 << '\n';
   // transform1
   std::cout << "translate pA to origin\n";
#endif
   auto const pA1 = pA0 - pA0;
   auto const pB1 = pB0 - pA0;
   auto const pC1 = pC0 - pA0;
#if defined SHOW_MATRIX_CALC
//#######################################################################
   auto mt = quan::fusion::make_translation_matrix(-A.centre);

   display(mt,"translation_matrix");



   auto pB1v = pB0v * mt   ;
//#############################################################
#endif

#if defined DEBUG_PRINT
   std::cout << "pA1 = " << pA1 << '\n';
   std::cout << "pB1 = " << pB1 << '\n';
#if defined SHOW_MATRIX_CALC
   display(pB1v, "pB1v = ");
#endif
   std::cout << "pC1 = " << pC1 << '\n';
#endif
   assert ( (pA1 == point{0_km,0_km,0_km}) );

   // transform 2 ---------------------------------
   assert( abs(pB1.x) > epsilon_km);
   auto const y_angle = quan::angle::deg(quan::atan2(pB1.z,pB1.x));
#if defined DEBUG_PRINT
   std::cout << "rotate around y-axis by " << y_angle << " so that pB.z == 0)\n";
#endif
   quan::three_d::y_rotation y_rotate{-y_angle};

#if defined SHOW_MATRIX_CALC
   //######################################################
   auto mry = quan::fusion::make_3d_y_rotation_matrix<quan::length::km>(-y_angle);
   auto pB2v = pB1v * mry;
   //######################################################
#endif
   auto const pA2 = y_rotate(pA1); 
   auto const pB2 = y_rotate(pB1); 

   auto const pC2 = y_rotate(pC1);

#if defined DEBUG_PRINT
   std::cout << "pA2 = " << pA2 << '\n';
   std::cout << "pB2 = " << pB2 << '\n';
#if defined SHOW_MATRIX_CALC
   display(pB2v, "pB2v = ");
#endif
   std::cout << "pC2 = " << pC2 << '\n';
#endif
   assert((pA2 == point{0_km,0_km,0_km}));
   assert(abs(pB2.z) < epsilon_km);

   // transform 3
   assert( abs(pB2.x) > epsilon_km);
   auto const z_angle = quan::angle::deg(quan::atan2(pB2.y,pB2.x));
#if defined DEBUG_PRINT
   std::cout << "rotate around z-axis by " << z_angle << " so that pB.y == 0\n";
#endif
   quan::three_d::z_rotation z_rotate{-z_angle};
#if defined SHOW_MATRIX_CALC
   //######################################################
   auto mrz = quan::fusion::make_3d_z_rotation_matrix<quan::length::km>(-z_angle);
   auto pB3v = pB2v * mrz;
   //######################################################
#endif

   auto const pA3 = z_rotate(pA2); 
   auto const pB3 = z_rotate(pB2); 
   auto const pC3 = z_rotate(pC2);

#if defined DEBUG_PRINT
   std::cout << "pA3 = " << pA3 << '\n';
   std::cout << "pB3 = " << pB3 << '\n';
#if defined SHOW_MATRIX_CALC
   display(pB3v, "pB3v = ");
#endif
   std::cout << "pC3 = " << pC3 << '\n';
#endif

   assert((pA3 == point{0_km,0_km,0_km}));
   assert(abs(pB3.z) < epsilon_km);
   assert(abs(pB3.y) < epsilon_km); 

   assert( abs(pC3.x) > epsilon_km);
   auto const x_angle = quan::angle::deg(quan::atan2(pC3.z,pC3.y));
#if defined DEBUG_PRINT
   std::cout << "rotate around x-axis by " << x_angle << " so that pC.z == 0\n";
#endif
   quan::three_d::x_rotation x_rotate{-x_angle};

#if defined SHOW_MATRIX_CALC
   //################################################################################

   auto mrx = quan::fusion::make_3d_x_rotation_matrix<quan::length::km>(-x_angle);
   display(mrx, "mrx = " ) ;

   auto pB4v = pB3v * mrx  ;

   auto mxtot = mt * mry * mrz * mrx;

   display(mxtot,"mxtot = ");

   auto pb4vv = pB0v * mxtot;
  //###################################################################################
#endif
   auto const pA4 = x_rotate(pA3); 
   auto const pB4 = x_rotate(pB3); 
   auto const pC4 = x_rotate(pC3);

#if defined DEBUG_PRINT
   std::cout << "pA4 = " << pA4 << '\n';
   std::cout << "pB4 = " << pB4 << '\n';
#if defined SHOW_MATRIX_CALC
   display(pB4v, "pB4v = ");
   display(pb4vv,"pb4vv ( calced from totala matrix) = ");
#endif
   std::cout << "pC4 = " << pC4 << '\n';
#endif

   assert((pA4 == point{0_km,0_km,0_km}));
   assert(abs(pB4.z) < epsilon_km);
   assert(abs(pB4.y) < epsilon_km); 
   assert(abs(pC4.z) < epsilon_km);

   point ip4;

   if ( !ll_trilaterate(sphere{pA4,A.radius},sphere{pB4,B.radius},sphere{pC4,C.radius},ip4)){
      return false;
   }

#if defined SHOW_MATRIX_CALC
   //#######################################################################
   auto mrx_dash = quan::fusion::make_3d_x_rotation_matrix<quan::length::km>(x_angle);
   auto mrz_dash = quan::fusion::make_3d_z_rotation_matrix<quan::length::km>(z_angle);
   auto mry_dash = quan::fusion::make_3d_y_rotation_matrix<quan::length::km>(y_angle);
   auto mt_dash = quan::fusion::make_translation_matrix(A.centre);

   auto mxtot_dash = mrx_dash * mrz_dash * mry_dash * mt_dash;

   auto ip4v = quan::fusion::make_column_matrix(ip4);

   auto ip0v = ip4v * mxtot_dash;
#if defined DEBUG_PRINT
   display(ip0v,"ip0v = ");
#endif

#endif
   //#######################################################################

   quan::three_d::x_rotation x_unrotate(x_angle);

   point const ip3 = x_unrotate(ip4);
#if defined DEBUG_PRINT
   point const pA3_dash = x_unrotate(pA4);
   point const pB3_dash = x_unrotate(pB4);
   point const pC3_dash = x_unrotate(pC4);
#endif
   quan::three_d::z_rotation z_unrotate(z_angle);
   point const ip2 = z_unrotate(ip3);
#if defined DEBUG_PRINT
   point pA2_dash = z_unrotate(pA3_dash);
   point pB2_dash = z_unrotate(pB3_dash);
   point pC2_dash = z_unrotate(pC3_dash);
#endif
   quan::three_d::y_rotation y_unrotate(y_angle);
   point ip1 = y_unrotate(ip2);
#if defined DEBUG_PRINT
   point pA1_dash = y_unrotate(pA2_dash);
   point pB1_dash = y_unrotate(pB2_dash);
   point pC1_dash = y_unrotate(pC2_dash);
#endif
   point ip0 = ip1 + pA0;
#if defined DEBUG_PRINT
   point pA0_dash = pA1_dash + pA0;
   point pB0_dash = pB1_dash + pA0;
   point pC0_dash = pC1_dash + pA0; 
   std::cout << "pA0_dash = " << pA0_dash << '\n';
   std::cout << "pB0_dash = " << pB0_dash << '\n';
   std::cout << "pC0_dash = " << pC0_dash << '\n';
#endif
   out = ip0;
   return true;
}

void output_scad_preamble(std::ostream & out)
{
   out << "// OpenScad script\n";
   out << "// https://www.openscad.org/\n\n";
   out << "// Generated by \"trilateration_transform.cpp\"\n";
   out << "// https://github.com/kwikius/Trilateration\n\n";
   out << "module show_sphere(pos,radius)\n";
   out << "{\n";
   out << "   translate(pos){\n";
   out << "      sphere(r = radius, $fn = 50);\n";
   out << "   }\n";
   out << "}\n\n";
// generated by "trilateration_transform.cpp"
}

int main()
{
  
   sphere A{{4_km, 5_km,6_km},7.5_km};;
   sphere B{{13_km, 4.5_km, 5.5_km},5.0_km};
   sphere C{{10_km,11_km,5.6_km},7.0_km};
#if defined DEBUG_PRINT
   std::cout << "A = " << A <<'\n';
   std::cout << "B = " << B <<'\n';
   std::cout << "C = " << C <<'\n';
#endif
   
   point intersection_point;
   if ( trilaterate(A,B,C,intersection_point)){
#if defined DEBUG_PRINT
     std::cout << "intersection point = " << intersection_point << "\n";

#endif
     {
        std::ofstream out("trilateration_transform.scad");

        output_scad_preamble(out);

        out << "color(\"blue\"){\n";
        out << "   // sphere A\n";
        out << "   show_sphere(" << A.centre/1_km << ", " << A.radius.numeric_value() << ");\n";
       out <<  "   // sphere B\n";
        out << "   show_sphere(" << C.centre/1_km << ", " << C.radius.numeric_value() << ");\n";
       out <<  "   // sphere C\n";
        out << "   show_sphere(" << B.centre/1_km << ", " << B.radius.numeric_value() << ");\n"; 
        out << "}\n\n";

        out << "color(\"yellow\"){\n";
        out << "   // sphere at intersection point\n";
        out << "   show_sphere(" << intersection_point / 1_km << ", 1);\n";
        out << "}\n\n";
     }

     std::cout << "result output to \"trilateration_transform.scad\"\n"; 
     std::cout << "...opening in openscad" << std::endl;

     return  system("openscad trilateration_transform.scad");

   }else{
      std::cout << "failed to trilaterate\n";
   }

   return 0;
}
