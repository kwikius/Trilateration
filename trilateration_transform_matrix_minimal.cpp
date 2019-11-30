
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

//#define SHOW_VECT_CALC
//#define SHOW_MATRIX_CALC

//#define USE_MATRIX_CALC
#define USE_VECT_CALC

#if defined (USE_VECT_CALC) && defined(USE_MATRIX_CALC)
#error choose calc
#endif

#if ! (defined (USE_VECT_CALC) || defined(USE_MATRIX_CALC))
#error need 1 or other calc
#endif

#if defined(SHOW_MATRIX_CALC) || defined (USE_MATRIX_CALC)
#define WANT_MATRIX_CALC
#endif

#if defined(SHOW_VECT_CALC) || defined (USE_VECT_CALC)
#define WANT_VECT_CALC
#endif



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

      auto const ex = unit_vector(B.centre);   // direction of B to origin
      auto const i = dot_product(ex,(C.centre));   
      auto const ey = unit_vector(C.centre - i * ex) ;
      auto const d = magnitude(B.centre);       // distance B to origin
      auto const j = dot_product(ey,C.centre);
      auto const x = (quan::pow<2>(A.radius) - quan::pow<2>(B.radius) + quan::pow<2>(d)) / ( 2 * d);

      if ( ( (d - A.radius) >= B.radius ) || ( B.radius >= (d + A.radius) ) ){
         // shouldnt get here as was checked in parent function
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
   //      std::cout << "ll intersection point = " << intersection_point << '\n';
   #endif
         return true;
      }else{
         std::cout << "z : no solution\n";
         return false;
      }
      return true;
   }

}

namespace quan{ namespace fun{
    
   template <typename M>
   inline
   typename quan::where_< 
      quan::meta::and_<
         quan::fun::is_fun_matrix<M>
         ,quan::meta::bool_< (quan::fusion::num_rows<M> == 1)>
         ,quan::meta::bool_<(quan::fusion::num_columns<M> == 4)>
      >,
      quan::three_d::vect<typename quan::fusion::matrix_at_t<0,0,M>::type>
   >::type
   to_vect3D( M const & in)
   {
       typedef typename quan::fusion::matrix_at_t<0,0,M>::type value_type;
       return quan::three_d::vect<value_type>{in. template at<0,0>(),in. template at<0,1>(),in. template at<0,2>()};
   }

}}

bool trilaterate_verify(sphere const& A, sphere const & B, sphere const & C)
{
   auto const distAB = magnitude(A.centre-B.centre);
   if (  distAB < epsilon_km ){
      std::cout << "A and B are coincident\n";
      return false;
   }
   if ( distAB >= (A.radius + B.radius) ){
      std::cout << "A and B dont intersect\n";
      return false;
   }
   auto const distBC = magnitude(B.centre-C.centre);
   if (  distBC < epsilon_km ){
      std::cout << "B and C are coincident\n";
      return false;
   }
   if ( distBC >= (B.radius + C.radius) ){
      std::cout << "B and C dont intersect\n";
      return false;
   }
   auto const distAC = magnitude(A.centre-C.centre);
   if ( distAC  < epsilon_km ){
      std::cout << "A and C are coincident\n";
      return false;
   }
   if ( distAC >= (A.radius + C.radius) ){
      std::cout << "A and C dont intersect\n";
      return false;
   }
   return true;
        
}
/*
  to align for matrix calc
  B1, C1  <- translate B,C by -A.centre
  mt 
  read y_angle  ( atan2(B1.z,B1.x))
  B2, C2 <- rotate B1,C1 around y by yangle
  my
  read z_angle (atan2(pB2.y, pB2.x)
  C3 <-- rotate C2 around z by z_angle
  mz
  read x_angle ( atan2(C3.z, C3.y)
   
  // that gets angles to rotate by
  // if want matrix then make matrix mt * my * mz * mx
  // apply to point
*/

bool trilaterate(sphere const& A, sphere const & B, sphere const & C,point & out)
{
   if (!trilaterate_verify(A,B,C)){
      return false;
   }

#if defined WANT_MATRIX_CALC
   auto const pA0v = quan::fusion::make_column_matrix(A.centre);
   auto const pB0v = quan::fusion::make_column_matrix(B.centre);
   auto const pC0v = quan::fusion::make_column_matrix(C.centre);
#endif

#if defined DEBUG_PRINT
   std::cout << "\ntranslate system so that A is at origin --------------\n\n";
#endif

#if defined WANT_VECT_CALC
   auto const pA_norm = A.centre - A.centre;
   assert ( (pA_norm == point{0_km,0_km,0_km}) );
   auto const pB1 = B.centre - A.centre;
   assert( abs(pB1.x) > epsilon_km);
   auto const pC1 = C.centre - A.centre;
#endif
#if defined WANT_MATRIX_CALC
   auto mt = quan::fusion::make_translation_matrix(-A.centre);
   auto const pAv_norm = pA0v * mt   ;
   auto const pB1v = pB0v * mt   ;
   auto const pC1v = pC0v * mt   ;
#endif

#if defined DEBUG_PRINT
   #if defined SHOW_VECT_CALC
      std::cout << "pA_norm = " << pA_norm << '\n';
   #endif
   #if defined SHOW_MATRIX_CALC
      display(pAv_norm, "pAv_norm = ");
   #endif
    #if defined SHOW_VECT_CALC
   std::cout << "pB1 = " << pB1 << '\n';
   #endif
   #if defined SHOW_MATRIX_CALC
      display(pB1v, "pB1v = ");
   #endif
   #if defined SHOW_VECT_CALC
       std::cout << "pC1 = " << pC1 << '\n';
   #endif
   #if defined SHOW_MATRIX_CALC
      display(pC1v, "pC1v = ");
   #endif
#endif   
   // transform 2 ---------------------------------
#if defined USE_MATRIX_CALC
   assert( (pB1v.at<0,3>()== 1) );
   auto const y_angle = quan::atan2(pB1v.at<0,2>(), pB1v.at<0,0>());
#else
   auto const y_angle = quan::atan2(pB1.z,pB1.x);
#endif
#if defined DEBUG_PRINT
   std::cout << "\nrotate around y-axis by " << quan::angle::deg{y_angle} << " so that pB.z == 0)\n\n";
#endif
#if defined WANT_MATRIX_CALC
   auto const mry = quan::fusion::make_3d_y_rotation_matrix<quan::length::km>(-y_angle);
#endif
#if defined WANT_VECT_CALC
   quan::three_d::y_rotation y_rotate{-y_angle};
#endif

#if defined WANT_MATRIX_CALC
   auto const pB2v = pB1v * mry;
   auto const pC2v = pC1v * mry;
#endif
#if defined WANT_VECT_CALC
   auto const pB2 = y_rotate(pB1); 
   auto const pC2 = y_rotate(pC1);
#endif
#if defined DEBUG_PRINT
   #if defined SHOW_VECT_CALC
      std::cout << "pB2 = " << pB2 << '\n';
   #endif
   #if defined SHOW_MATRIX_CALC
      display(pB2v, "pB2v = ");
   #endif
   #if defined SHOW_VECT_CALC
       std::cout << "pC2 = " << pC2 << '\n';
   #endif
   #if defined SHOW_MATRIX_CALC
      display(pC2v, "pC2v = ");
   #endif
#endif
  
#if defined USE_MATRIX_CALC
   assert( (pB2v.at<0,3>()== 1) );
   auto const z_angle = quan::atan2(pB2v.at<0,1>(),pB2v.at<0,0>());
#else
   auto const z_angle = quan::atan2(pB2.y,pB2.x);
#endif

#if defined DEBUG_PRINT
   std::cout << "\nrotate around z-axis by " << quan::angle::deg{z_angle} << " so that pB.y == 0\n\n";
#endif

#if defined WANT_MATRIX_CALC
   auto mrz = quan::fusion::make_3d_z_rotation_matrix<quan::length::km>(-z_angle);
   auto pBv_norm = pB2v * mrz;
   auto pC3v = pC2v * mrz;
#endif
#if defined WANT_VECT_CALC
   quan::three_d::z_rotation z_rotate{-z_angle};
   auto const pB_norm = z_rotate(pB2); 
   assert(abs(pB_norm.z) < epsilon_km);
   assert(abs(pB_norm.y) < epsilon_km);
 
   auto const pC3 = z_rotate(pC2);
   assert( abs(pC3.x) > epsilon_km);
#endif

#if defined DEBUG_PRINT
   #if defined SHOW_VECT_CALC
      std::cout << "pB_norm = " << pB_norm << '\n';
   #endif
   #if defined SHOW_MATRIX_CALC
      display(pBv_norm, "pBv_norm = ");
   #endif
   #if defined SHOW_VECT_CALC
      std::cout << "pC3 = " << pC3 << '\n';
   #endif
   #if defined SHOW_MATRIX_CALC
      display(pC3v, "pC3v = ");
   #endif
#endif
 

#if defined WANT_MATRIX_CALC
   assert( (pC3v.at<0,3>()== 1) );
   auto const x_angle = quan::atan2(pC3v.at<0,2>(),pC3v.at<0,1>());
#else
   auto const x_angle = quan::atan2(pC3.z,pC3.y);
#endif

#if defined DEBUG_PRINT
   std::cout << "\nrotate around x-axis by " << quan::angle::deg{x_angle} << " so that pC.z == 0\n\n";
#endif

#if defined WANT_MATRIX_CALC
   auto mrx = quan::fusion::make_3d_x_rotation_matrix<quan::length::km>(-x_angle);
   #if defined SHOW_MATRIX_CALC
      display(mrx, "mrx = " ) ;
   #endif
#endif
#if defined WANT_VECT_CALC
   quan::three_d::x_rotation x_rotate{-x_angle};
   auto const pC_norm = x_rotate(pC3);
   assert(abs(pC_norm.z) < epsilon_km);
   #if defined SHOW_VECT_CALC
      std::cout << "pC_norm = " << pC_norm << '\n';
   #endif
#endif
#if defined WANT_MATRIX_CALC
   auto const pCv_norm = pC3v * mrx;
#endif

   point ip_norm;
#if defined WANT_MATRIX_CALC
   if ( !ll_trilaterate(
             sphere{to_vect3D(pAv_norm),A.radius}
            ,sphere{to_vect3D(pBv_norm),B.radius}
            ,sphere{to_vect3D(pCv_norm),C.radius}
            ,ip_norm
         )
   ){
      return false;
   }
#else
   if ( !ll_trilaterate(sphere{pA_norm,A.radius},sphere{pB_norm,B.radius},sphere{pC_norm,C.radius},ip_norm)){
      return false;
   }
#endif

#if defined USE_MATRIX_CALC
   auto mrx_dash = quan::fusion::make_3d_x_rotation_matrix<quan::length::km>(x_angle);
   auto mrz_dash = quan::fusion::make_3d_z_rotation_matrix<quan::length::km>(z_angle);
   auto mry_dash = quan::fusion::make_3d_y_rotation_matrix<quan::length::km>(y_angle);
   auto mt_dash = quan::fusion::make_translation_matrix(A.centre);

   auto mxtot_dash = mrx_dash * mrz_dash * mry_dash * mt_dash;

   auto ipv_norm = quan::fusion::make_column_matrix(ip_norm);
   auto ip0v = ipv_norm * mxtot_dash;
   out = to_vect3D(ip0v);
#else
   quan::three_d::x_rotation x_unrotate(x_angle);
   point const ip3 = x_unrotate(ip_norm);
   quan::three_d::z_rotation z_unrotate(z_angle);
   point const ip2 = z_unrotate(ip3);
   quan::three_d::y_rotation y_unrotate(y_angle);
   point ip1 = y_unrotate(ip2);
   point ip0 = ip1 + A.centre;
   out = ip0;
#endif
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
     std::cout << "\nintersection point = " << intersection_point << "\n\n";

#endif
     {
        std::ofstream out("trilateration_transform.scad");

        output_scad_preamble(out);

        out << "color(\"blue\"){\n";
        out << "   // sphere A\n";
        out << "   show_sphere(" << A.centre/1_km << ", " << A.radius.numeric_value() << ");\n";
        out << "   // sphere B\n";
        out << "   show_sphere(" << C.centre/1_km << ", " << C.radius.numeric_value() << ");\n";
        out << "   // sphere C\n";
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
