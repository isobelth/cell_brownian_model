#include "colors.inc"
camera {
  sky <0,0,1>           
  direction <-1,0,0>      
  right <-4/3,0,0>      
  location <0,-5,10> 
  look_at <0,0,0>     
  angle 15      
}
global_settings { ambient_light White }
light_source {
  <10,-10,20>   
  color White*2 
}
background { color White }
