
/*
 * Implementation of Stochastic Yan Network using Gillespie framework
 */

#include <iomanip>

#include "basic_gillespie.h"

class Degradation : public Reaction {
public:
  Degradation(std::string sp,double rate) : species(sp), k(rate)
  {
    // initializing base class fields
    name = sp + "_deg";
    parameters = Params({{"k_"+sp,rate}});
  }
  double propensity(State &s,double t) {return k * s[species];}
  void perform(State &s) {--s[species];}
private:
  std::string species;
  double k;
};

class ConstInduction : public Reaction {
public:
  ConstInduction(std::string sp,double rate) : species(sp), alpha(rate)
  {
    // initializing base class fields
    name = sp + "_ind";
    parameters = Params({{"alpha_"+sp,rate}});
  }
  double propensity(State &s,double t) {return alpha;}
  void perform(State &s) {++s[species];}
private:
  std::string species;
  double alpha;
};

class TimeInduction : public Reaction {
public:
  TimeInduction(std::string sp,std::function<double(double)> ind_fun) : species(sp), alpha(ind_fun)
  {
    // initializing base class fields
    name = sp + "_ind";
    parameters = Params({});
  }
  double propensity(State &s,double t) {return alpha(t);}
  void perform(State &s) {++s[species];}
private:
  std::string species;
  double time;
  std::function<double(double)> alpha;
};

class PRTInduction : public Reaction {
public:
  PRTInduction(std::string sp,std::string msp,double rate) : species(sp), mspecies(msp),k(rate)
  {
    name = sp + "_ind";
    parameters = Params({{"k_"+sp,rate}});
  }
  double propensity(State &s,double t) {return k * s[mspecies];}
  void perform(State &s) {++s[species];}
private:
  std::string species,mspecies;
  double k;
};

class HillInduction : public Reaction{
public:
  HillInduction(std::string sp,std::string hsp,double max_rate,double k,double n) : HillInduction(sp,std::vector<std::string>{hsp},max_rate,k,n) {}
  // for activation by a set of species
  HillInduction(std::string sp,std::vector<std::string> hsps,double max_rate,double k,double n) : species(sp), activators(hsps), alpha(max_rate), K(k), N(n)
  {
    name = sp + "_ind";
    parameters = Params({{"alpha_" + sp,alpha},{"K",K},{"N",N}});
  }  
  double propensity(State &s,double t)
  {
    double total_activator = get_activators(s);
    return alpha * pow(total_activator,N) / (pow(K,N) + pow(total_activator,N));
  }
  void perform(State &s) {++s[species];}
private:
  double get_activators(State &s)
  {
    double total = 0.0;
    for (auto a : activators) {total += s[a];}
    return total;
  }
  std::string species;
  std::vector<std::string> activators;
  double alpha,K,N;
};

class YPSiteInduction : public Reaction {
public:
  YPSiteInduction(std::string sp,std::vector<double> &rates, std::vector< std::pair<std::string,int> > &mults, std::vector<double> &dGs,double vol,double t) : species(sp), alphas(rates), multiplicities(mults), deltaGs(dGs), V(vol), time(t)
  {
    name = sp + "_ind";
    parameters = Params();
  }
  double propensity(State &s,double t)
  {
    std::vector<double> bweights(alphas.size(),0.0);
    for (size_t i = 0; i < alphas.size(); ++i) {
      bweights[i] = calc_weight(multiplicities[i],deltaGs[i],s);
    }
    double btotal = std::accumulate(bweights.begin(),bweights.end(),0.0);
    double eff_rate = 0.0;
    for (size_t i = 0; i < alphas.size(); ++i) {
      eff_rate += alphas[i] * bweights[i] / btotal;
    }
    return eff_rate * (t > time);
  }
  void perform(State &s) {++s[species];}
private:
  double calc_conc(int copy_number) {return copy_number / (6.022e23 * V);}
  double calc_weight(std::pair<std::string,int> m,double dG,State &s)
  {
    double KbT = 0.593;
    if (m.first == "")
      return exp(dG/KbT);
    return pow(calc_conc(s[m.first]),m.second) * exp(dG/KbT);
  }
  std::string species;
  std::vector<double> alphas;
  std::vector< std::pair<std::string,int> > multiplicities;
  std::vector<double> deltaGs;
  double V;
  double time;
};

class TimePhosphorylation : public Reaction {
public:
  TimePhosphorylation(std::string sp,std::string psp,std::function<double(double)> efun,double max_rate,double k) : species(sp), pspecies(psp), E(efun), alpha(max_rate), K(k)
  {
    name = sp + "_phospho";
    parameters = Params({{"alpha_" + sp,alpha},{"K",K}});
  }
  double propensity(State &s,double t)
  {
    return alpha * E(t) * s[species] / (K + s[species]);
  }
  void perform(State &s) {++s[pspecies];--s[species];}
private:
  std::string species,pspecies;
  double alpha,K,start_time,end_time;
  std::function<double(double)> E;
};

class ComplexFormation : public Reaction {
public:
  ComplexFormation(std::string sp1,std::string sp2,std::string cmp,double rate) : species1(sp1), species2(sp2), complex(cmp), k(rate)
  {
    name = complex +  "_form";
    parameters = Params{{"k_" + complex,rate}};
  }
  double propensity(State &s,double t) {return k * s[species1] * s[species2];}
  void perform(State &s) {++s[complex];--s[species1];--s[species2];}
private:
  std::string species1,species2,complex;
  double k;
};

class ComplexDissociation : public Reaction {
public:
  ComplexDissociation(std::string cmp,std::string sp1,std::string sp2,double rate) : species1(sp1), species2(sp2), complex(cmp), k(rate)
  {
    name = complex + "_diss";
    parameters = Params{{"k_" + complex,rate}};
  }
  double propensity(State &s,double t) {return k * s[complex];}
  void perform(State &s) {--s[complex];++s[species1];++s[species2];}
private:
  std::string species1,species2,complex;
  double k;
};

/*
 * Add all reactions here in Yan network
 */
void yan_network_init(BasicGillespie &BG) {
  // General Constants
  double S2_NUC = 78e-15; // 78 um^3 volume in litres
  
  // Degradation
  // mRNA
  double DEF_RNA_DEG = 1.0/600; // ~10min mean lifetime default mRNA degradation
  BG.add_reaction(std::make_shared<Degradation>("mY",DEF_RNA_DEG));
  BG.add_reaction(std::make_shared<Degradation>("mR7",DEF_RNA_DEG));
  BG.add_reaction(std::make_shared<Degradation>("mR7_mY",DEF_RNA_DEG));
  
  // Protein
  double DEF_PRT_DEG = 1.0/(5.0*3600); // ~1hr mean lifetime default protein degradation
  double YP_PRT_DEG = 1.0/300; // ~5 min lifetime for phosphorylated yan
  BG.add_reaction(std::make_shared<Degradation>("Y",DEF_PRT_DEG));
  BG.add_reaction(std::make_shared<Degradation>("YP",YP_PRT_DEG));
  BG.add_reaction(std::make_shared<Degradation>("P2",DEF_PRT_DEG));
  BG.add_reaction(std::make_shared<Degradation>("P2P",DEF_PRT_DEG));
  BG.add_reaction(std::make_shared<Degradation>("P1",DEF_PRT_DEG));
  BG.add_reaction(std::make_shared<Degradation>("M",DEF_PRT_DEG));
  double SUH_PRT_DEG = 1.0/3600; // ~1hr, faster response
  BG.add_reaction(std::make_shared<Degradation>("S",SUH_PRT_DEG)); // Su(H)
  BG.add_reaction(std::make_shared<Degradation>("Y_Y",DEF_PRT_DEG));
  BG.add_reaction(std::make_shared<Degradation>("M_P2P",DEF_PRT_DEG));
  BG.add_reaction(std::make_shared<Degradation>("M_Y",DEF_PRT_DEG));
  BG.add_reaction(std::make_shared<Degradation>("M_YP",DEF_PRT_DEG));

  // Some energy constants for promoter models (Mae,miR7)
  double dGE = 10.0; // 10.0 original
  double dGNS = 5.8;
  double dGSAM = 7.0;
  double dGP = 12.4; // 8.2 original --> 12.4 high
  // for Yan binding site
  double dGEmY = 8.0; // weaker yan ets binding
  double dGPmY = 14.0; // strongest pnt ets binding
  
  // RNA Induction
  double SUH_PRT_COPY = 1;
  double DEF_RNA_COPY = 20;
  double mY_start = 0*3600.0;
  // mY RNA controlled explicitly in time
  // std::function<double(double)> YAN_IND = [&](double t) {return DEF_RNA_DEG*DEF_RNA_COPY *(t > mY_start);};
  // BG.add_reaction(std::make_shared<TimeInduction>("mY",YAN_IND));
  // mY controlled by Su(H)
  // BG.add_reaction(std::make_shared<PRTInduction>("mY","S",DEF_RNA_DEG*DEF_RNA_COPY/SUH_PRT_COPY));
  // Yan promoter binding site model
  double YAN_LOW_COPY = 50; // baseline expression copy number
  double YAN_PRT_COPY = 500; // high state Yan copy number
  double mY_max = DEF_RNA_DEG*DEF_RNA_COPY * 1.5;
  double mY_low = DEF_RNA_DEG*DEF_RNA_COPY*YAN_LOW_COPY/YAN_PRT_COPY;
  // 8 states now with Su(H) binding
  std::vector<double> mY_rates{mY_low,mY_low,mY_low,0,0,mY_max,mY_max,mY_max};
  std::vector< std::pair<std::string,int> > mY_mults =
    {std::pair<std::string,int> {"",0},
     std::pair<std::string,int> {"Y",1},
     std::pair<std::string,int> {"Y",1},
     std::pair<std::string,int> {"Y",2},
     std::pair<std::string,int> {"Y_Y",1},
     std::pair<std::string,int> {"P1",1},
     std::pair<std::string,int> {"P2P",1},
     std::pair<std::string,int> {"S",1}};
  std::vector<double> mY_dGs {0,dGEmY,dGNS,dGEmY+dGNS+dGSAM,dGEmY+dGNS,dGPmY,dGPmY,dGPmY};
  BG.add_reaction(std::make_shared<YPSiteInduction>("mY",mY_rates,mY_mults,mY_dGs,S2_NUC,mY_start));
  BG.add_pause(mY_start);

  // Protein Induction
  double DEF_PRT_COPY = 500;
  // Yan Induction
  BG.add_reaction(std::make_shared<PRTInduction>("Y","mY",DEF_PRT_DEG*YAN_PRT_COPY/DEF_RNA_COPY));
  // PntP2 Induction
  BG.add_reaction(std::make_shared<ConstInduction>("P2",DEF_PRT_DEG*DEF_PRT_COPY));
  // PntP1 Hill Induction
  double P1_rate = DEF_PRT_DEG*DEF_PRT_COPY*1;
  double P1_K = DEF_PRT_COPY/2.0;
  BG.add_reaction(std::make_shared<HillInduction>("P1","P2P",P1_rate,P1_K,1.0));
  // Alternative PntP1 Hill Induction with positive feedback
  // BG.add_reaction(std::make_shared<HillInduction>("P1",std::vector<std::string>{"P2P","P1"},DEF_PRT_DEG*DEF_PRT_COPY,DEF_PRT_COPY/3.0,2.0));
  // Su(H) const induction
  double SUH_start = 10*3600.0;
  std::function<double(double)> SUH_IND = [&](double t) {return SUH_PRT_DEG*SUH_PRT_COPY *(t > SUH_start);};
  BG.add_reaction(std::make_shared<TimeInduction>("S",SUH_IND));
  
  // Mae Binding Site Induction
  double MAE_PRT_COPY = 500;
  double M_max = DEF_PRT_DEG*MAE_PRT_COPY;
  double M_low = 0;
  double M_start = 15 * 3600.0; // same as erk_start
  std::vector<double> M_rates{M_low,M_low,M_low,0,0,M_max,M_max};
  std::vector< std::pair<std::string,int> > M_mults =
    {std::pair<std::string,int> {"",0},
     std::pair<std::string,int> {"Y",1},
     std::pair<std::string,int> {"Y",1},
     std::pair<std::string,int> {"Y",2},
     std::pair<std::string,int> {"Y_Y",1},
     std::pair<std::string,int> {"P1",1},
     std::pair<std::string,int> {"P2P",1}};
  std::vector<double> M_dGs{0,dGE,dGNS,dGE+dGNS+dGSAM,dGE+dGNS,dGP,dGP};
  BG.add_reaction(std::make_shared<YPSiteInduction>("M",M_rates,M_mults,M_dGs,S2_NUC,M_start));
  // mR7 Binding Site Induction -- using parameters as Mae induction
  double mR7_RNA_COPY = 100; // increased suppression of Yan
  double mR7_max = DEF_RNA_DEG*mR7_RNA_COPY;
  double mR7_low = 0;
  double mR7_start = M_start;
  std::vector<double> mR7_rates{mR7_low,mR7_low,mR7_low,0,0,mR7_max,mR7_max};
  std::vector< std::pair<std::string,int> > mR7_mults = M_mults;
  std::vector<double> mR7_dGs = M_dGs;
  BG.add_reaction(std::make_shared<YPSiteInduction>("mR7",mR7_rates,mR7_mults,mR7_dGs,S2_NUC,mR7_start));
  
  // Phophorylation
  double DEF_PRT_PHOSPHO = 1.0/60; // ~1 min time to phosphorylate per ERK
  double P1_P2P_target = 20;
  double ERK_HIGH = 4.0;
  double ERK_LOW = P1_P2P_target / ((DEF_PRT_PHOSPHO/DEF_PRT_DEG)*(1 + P1_rate/(DEF_PRT_DEG*P1_K))); // 10.0 high?
  std::cout << ERK_LOW << std::endl;
  double ERK_MIN = 0.0;
  double erk_start = 10 * 3600.0;
  double erk_pulse_start = 20 * 3600;
  double erk_pulse_end = 25 * 3600.0;
  double erk_end = 20 * 3600.0;
  std::function<double(double)> ERK_FUN = [=](double t)
    {
      if (t < erk_start)
	return 0.0;
      else if (t >= erk_start && t < erk_pulse_end)
	return ERK_LOW;
      else if (t >= erk_pulse_end && t < erk_end) {
	double a = (ERK_LOW - ERK_MIN) / (erk_pulse_end - erk_end);
	double b = ERK_LOW - a * erk_pulse_end;
	return a * t + b;
      }
      else
	return 0.0;
    };
  // ERK function for pulse
  std::function<double(double)> ERK_PULSE = [=](double t)
    {
      if (t >= erk_pulse_start && t < erk_pulse_end)
	return ERK_HIGH;
      if (t > erk_pulse_end)
	return 0.5; // 0.5 regular extended
      else
	return ERK_FUN(t);
    };
  
  // PntP2 phosphorylation
  double ERK_SCALE = 1.0; // low scale 0.3 high scale 2.0
  std::function<double(double)> ERK_COURSE = [=](double t){return ERK_FUN(t)*ERK_SCALE;};
  BG.add_reaction(std::make_shared<TimePhosphorylation>("P2","P2P",ERK_COURSE,DEF_PRT_PHOSPHO,DEF_PRT_COPY/100));
  BG.add_pause(erk_start);
  // Mae:Yan phosphorylation
  BG.add_reaction(std::make_shared<TimePhosphorylation>("M_Y","M_YP",ERK_COURSE,DEF_PRT_PHOSPHO,MAE_PRT_COPY/3));
  
  // Interactions
  double DEF_CMP_DISS = 1.0/300; // ~5 min complex lifetime
  // Yan dimer formation
  double KdYY = 330000.0; // 7uM
  BG.add_reaction(std::make_shared<ComplexFormation>("Y","Y","Y_Y",DEF_CMP_DISS/KdYY));
  BG.add_reaction(std::make_shared<ComplexDissociation>("Y_Y","Y","Y",DEF_CMP_DISS));
  // Mae:PntP2P complex
  double KdMP2P = 470; // 10nM
  BG.add_reaction(std::make_shared<ComplexFormation>("M","P2P","M_P2P",DEF_CMP_DISS/KdMP2P));
  BG.add_reaction(std::make_shared<ComplexDissociation>("M_P2P","M","P2P",DEF_CMP_DISS));
  // Mae:Yan complex
  double KdMY = 470; // 10nM
  BG.add_reaction(std::make_shared<ComplexFormation>("M","Y","M_Y",DEF_CMP_DISS/KdMY));
  BG.add_reaction(std::make_shared<ComplexDissociation>("M_Y","M","Y",DEF_CMP_DISS));
  // Mae:YanP complex
  // more rapid dissociation
  double MY_CMP_DISS = 1.0/30; // ~30 sec complex lifetime
  double KdMYP = 4700; // 100nM
  BG.add_reaction(std::make_shared<ComplexFormation>("M","YP","M_YP",MY_CMP_DISS/KdMYP));
  BG.add_reaction(std::make_shared<ComplexDissociation>("M_YP","M","YP",MY_CMP_DISS));
  // mR7:mY complex
  double KdmR7mY = 10;
  BG.add_reaction(std::make_shared<ComplexFormation>("mR7","mY","mR7_mY",DEF_CMP_DISS/KdmR7mY));
  BG.add_reaction(std::make_shared<ComplexDissociation>("mR7_mY","mR7","mY",DEF_CMP_DISS));
}

int main(int argc,char const ** argv) {
  // Yan network initial conditions
  State IC{{"mY",2},{"mR7",0},{"mR7_mY",0},{"Y",50},{"YP",0},{"P2",500},{"P2P",0},{"P1",0},{"M",0},{"S",0},{"Y_Y",0},{"M_P2P",0},{"M_Y",0},{"M_YP",0}};

  // recording times
  std::vector<double> times;
  double dmax = 60.0 * 3600.0;
  for (double d = 0.0; d < dmax  || std::abs(d - dmax) < 1e-6; d += 60.0) {
    times.push_back(d);
  }

  std::string results_file("results/yan_network.txt");
  
  for (int i = 0; i < 3; ++i){
    BasicGillespie BG(IC,times);
    yan_network_init(BG);
    std::cout << "running " << i << std::endl;
    BG.run(dmax+1.0,1e6);
    // BG.print();
    if (i == 0)
      BG.write(results_file);
    else
      BG.write(results_file,false,true);
  }

  
}
