
#include "basic_gillespie.h"

int main(int argc, char const ** argv) {

  Reaction auto_x([](State &s,double t){return s["x"]*s["y"]*(t > 0.5);},
		  [](State &s){++s["x"];},
		  "x_auto",Params {{"k",1.0}});

  Reaction auto_y([](State &s,double t){return 1.0*(t > 0.5);},
		  [](State &s){++s["y"];},
		  "y_auto",Params {{"ky",1.0}});
  
  State IC{{"x",1},{"y",1}};
  std::vector<double> times {0.0,0.5,1.0,1.5,2.0,2.5};
  for (int i = 0; i < 3; ++i){
    BasicGillespie BG(IC,times);
    BG.add_reaction(std::make_shared<Reaction>(auto_x));
    BG.add_reaction(std::make_shared<Reaction>(auto_y));
    BG.add_pause(0.5);
    BG.run(3.0,1000);
    BG.print();
    if (i == 0)
      BG.write("results/testing.txt");
    else
      BG.write("results/testing.txt",false,true);
  }
  return 0;
}
